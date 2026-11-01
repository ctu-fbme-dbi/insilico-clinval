"""Synthetic novel loss-of-function variants for a test suite.

Two deterministic steps, VEP in between:

    python scripts/generate_synthetic.py generate \\
        --config configs/synthetic/exp_synthetic_v5q.yaml --out DIR
    # annotate DIR/candidates.vcf with VEP (pinned container, MANE Select, HGVS), e.g.
    #   docker run --rm -v $PWD:$PWD -w $PWD ensemblorg/ensembl-vep:release_115.0 vep \\
    #     --offline --cache --dir_cache data/vep/cache --assembly GRCh38 \\
    #     --fasta <reference FASTA> --mane_select --hgvs --vcf --force_overwrite \\
    #     --input_file DIR/candidates.vcf --output_file DIR/candidates.vep.vcf
    python scripts/generate_synthetic.py verify --config ... --out DIR \\
        --vep DIR/candidates.vep.vcf

`generate` builds nonsense, canonical splice-site and frameshift candidates on MANE Select
transcripts of the configured genes, keeps those predicted to trigger NMD and absent from ClinVar,
eRepo and (optionally) gnomAD, and ranks them by `derive_int(seed, "synthetic", key)`.
`verify` keeps candidates whose VEP consequence on the MANE Select transcript matches the intended
type, takes the target number per type and bin in rank order, assigns them to the positive samples
keeping `min_distance_bp` within each sample, and writes the suite from a template.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pysam
import yaml

from insilico_clinval.genome.alleles import left_normalize, ref_span0
from insilico_clinval.genome.bed import IntervalSet
from insilico_clinval.genome.reference import Reference
from insilico_clinval.rng import derive_int

STOPS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASES = "ACGT"
EXPECTED_CONSEQUENCE = {
    "nonsense": {"stop_gained"},
    "splice_donor": {"splice_donor_variant"},
    "splice_acceptor": {"splice_acceptor_variant"},
    "frameshift": {"frameshift_variant"},
}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# -- configuration ------------------------------------------------------------------------------


@dataclass
class Config:
    path: Path
    seed: int
    reference: Path
    mane_gtf: Path
    mane_release: str
    clingen_dosage: Path
    genes: list[str]
    exclude_genes: list[str]
    clinvar_vcf: Path
    erepo_tsv: Path | None
    erepo_release: str
    gnomad_vcf: Path | None
    avoid_strata: list[Path]
    avoid_pad: int
    cds_start_exclude: int
    nmd_margin: int
    indel_exon_margin: int
    min_distance_bp: int
    oversample: float
    max_per_gene: int
    targets: dict[str, int]
    frameshift_bins: dict[str, list[int]]
    frameshift_kinds: list[str]
    deletion_exclusion: str
    samples: list[str]
    template: Path
    suite_out: Path

    @staticmethod
    def load(path: Path) -> Config:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        base = path.parent

        def p(value: str | None) -> Path | None:
            return None if value is None else (base / value).resolve()

        def req(value: str) -> Path:
            out = p(value)
            assert out is not None
            return out

        return Config(
            path=path,
            seed=int(data["seed"]),
            reference=req(data["reference"]),
            mane_gtf=req(data["mane"]["gtf"]),
            mane_release=str(data["mane"]["release"]),
            clingen_dosage=req(data["clingen_dosage"]),
            genes=list(data["genes"]),
            exclude_genes=list(data.get("exclude_genes", [])),
            clinvar_vcf=req(data["clinvar_vcf"]),
            erepo_tsv=p(data.get("erepo", {}).get("file")),
            erepo_release=str(data.get("erepo", {}).get("release", "")),
            gnomad_vcf=p(data.get("gnomad_vcf")),
            avoid_strata=[req(s) for s in data.get("avoid_strata", [])],
            avoid_pad=int(data.get("avoid_pad", 5)),
            cds_start_exclude=int(data.get("cds_start_exclude", 100)),
            nmd_margin=int(data.get("nmd_margin", 50)),
            indel_exon_margin=int(data.get("indel_exon_margin", 3)),
            min_distance_bp=int(data["min_distance_bp"]),
            oversample=float(data.get("oversample", 1.5)),
            max_per_gene=int(data.get("max_per_gene", 0)),
            targets={str(k): int(v) for k, v in data["targets"].items()},
            frameshift_bins={
                str(k): [int(x) for x in v] for k, v in data["frameshift_bins"].items()
            },
            frameshift_kinds=list(data.get("frameshift_kinds", ["del", "ins"])),
            deletion_exclusion=str(data.get("deletion_exclusion", "span")),
            samples=list(data["samples"]),
            template=req(data["template"]),
            suite_out=req(data["suite_out"]),
        )

    def __post_init__(self) -> None:
        if not self.frameshift_kinds or set(self.frameshift_kinds) - {"del", "ins"}:
            raise ValueError(
                f"frameshift_kinds must be a subset of [del, ins]: {self.frameshift_kinds}"
            )
        if self.deletion_exclusion not in ("span", "breakpoints"):
            raise ValueError(
                f"deletion_exclusion must be span or breakpoints: {self.deletion_exclusion}"
            )


# -- gene models --------------------------------------------------------------------------------


@dataclass
class Transcript:
    gene: str
    transcript_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)  # 1-based closed, transcript order
    cds: list[tuple[int, int]] = field(default_factory=list)  # CDS without the stop codon

    def ordered(self) -> None:
        rev = self.strand == "-"
        self.exons.sort(key=lambda e: e[0], reverse=rev)
        self.cds.sort(key=lambda e: e[0], reverse=rev)

    def exon_of(self, pos: int) -> int:
        """0-based exon index (transcript order) containing genomic ``pos``; -1 if intronic."""
        for i, (s, e) in enumerate(self.exons):
            if s <= pos <= e:
                return i
        return -1

    def tx_coord(self, pos: int) -> int:
        """1-based mRNA coordinate of an exonic genomic position."""
        n = 0
        for s, e in self.exons:
            if s <= pos <= e:
                return n + (pos - s + 1 if self.strand == "+" else e - pos + 1)
            n += e - s + 1
        raise ValueError(f"{self.transcript_id}: {pos} is not exonic")

    def last_junction(self) -> int:
        """mRNA coordinate of the last base before the last exon-exon junction."""
        return sum(e - s + 1 for s, e in self.exons[:-1])

    def cds_positions(self) -> list[int]:
        """Genomic positions of the CDS bases in transcript order."""
        out: list[int] = []
        for s, e in self.cds:
            out.extend(range(s, e + 1) if self.strand == "+" else range(e, s - 1, -1))
        return out


def _attrs(field9: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for item in field9.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, value = item.partition(" ")
        out.setdefault(key, []).append(value.strip().strip('"'))
    return out


def read_mane(gtf: Path, genes: set[str]) -> dict[str, Transcript]:
    """MANE Select transcripts of the given genes (exons and CDS) from a MANE GTF."""
    found: dict[str, Transcript] = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("exon", "CDS"):
                continue
            a = _attrs(f[8])
            gene = a.get("gene_name", [""])[0]
            if gene not in genes or "MANE_Select" not in a.get("tag", []):
                continue
            t = found.setdefault(gene, Transcript(gene, a["transcript_id"][0], f[0], f[6]))
            (t.exons if f[2] == "exon" else t.cds).append((int(f[3]), int(f[4])))
    for t in found.values():
        t.ordered()
    return found


def haploinsufficient(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 4 and f[4] == "3":
                out.add(f[0])
    return out


# -- candidates ---------------------------------------------------------------------------------


@dataclass(frozen=True)
class Candidate:
    kind: str  # nonsense | splice_donor | splice_acceptor | frameshift
    bin: str  # "-" or the frameshift bin
    gene: str
    transcript: str
    chrom: str
    pos: int
    ref: str
    alt: str
    exon: int  # 1-based exon number in transcript order
    exons: int
    note: str

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def span(self) -> tuple[int, int]:
        return ref_span0(self.pos, self.ref)


class Exclusions:
    """Position-level exclusion: ClinVar, eRepo and gnomAD records overlapping a REF span.

    With ``deletion_exclusion: breakpoints`` a deletion is excluded only by records overlapping
    its breakpoints (the anchor and first deleted base, the last deleted base and the base after);
    an identical or larger record always does. In genes as densely covered by ClinVar as the ACMG
    SF genes, almost every deleted span of more than 2 bp holds some record, so the span rule
    leaves practically no deletions.
    """

    def __init__(self, cfg: Config, reference: Reference) -> None:
        self.reference = reference
        self.clinvar = pysam.VariantFile(str(cfg.clinvar_vcf))
        self.gnomad = pysam.VariantFile(str(cfg.gnomad_vcf)) if cfg.gnomad_vcf else None
        self.erepo: dict[str, list[tuple[int, int]]] = {}
        if cfg.erepo_tsv is not None:
            from insilico_clinval.variants.erepo import load_erepo

            loaded = load_erepo(
                cfg.erepo_tsv,
                release=cfg.erepo_release,
                reference=reference,
                held_out_reference=None,
            )
            for a in loaded.alleles:
                self.erepo.setdefault(a.chrom, []).append(ref_span0(a.pos, a.ref))
        self.strata = [IntervalSet(p) for p in cfg.avoid_strata]
        self.pad = cfg.avoid_pad
        self.breakpoints = cfg.deletion_exclusion == "breakpoints"
        self.counts: Counter[str] = Counter()

    @staticmethod
    def _contig(vcf: Any, chrom: str) -> str | None:
        names = set(vcf.header.contigs)
        for name in (chrom, chrom.removeprefix("chr"), "MT" if chrom == "chrM" else chrom):
            if name in names:
                return name
        return None

    def _vcf_hit(self, vcf: Any, chrom: str, s0: int, e0: int) -> bool:
        name = self._contig(vcf, chrom)
        if name is None:
            return False
        for rec in vcf.fetch(name, max(0, s0 - 1), e0 + 1):
            rs, re_ = ref_span0(rec.pos, str(rec.ref))
            if rs < e0 and s0 < re_:
                return True
        return False

    def reason(self, c: Candidate) -> str | None:
        s0, e0 = c.span
        if any(st.overlaps(c.chrom, max(0, s0 - self.pad), e0 + self.pad) for st in self.strata):
            return "repeat_context"
        windows = [(s0, e0)]
        if self.breakpoints and len(c.ref) > len(c.alt):
            windows = [(s0, s0 + 2), (e0 - 1, e0 + 1)]
        for ws, we in windows:
            if self._vcf_hit(self.clinvar, c.chrom, ws, we):
                return "in_clinvar"
        for ws, we in windows:
            if any(a < we and ws < b for a, b in self.erepo.get(c.chrom, [])):
                return "in_erepo"
        if self.gnomad is not None:
            for ws, we in windows:
                if self._vcf_hit(self.gnomad, c.chrom, ws, we):
                    return "in_gnomad"
        return None


class GeneSequence:
    """Reference sequence of one gene region held in memory (a ``SequenceSource``)."""

    def __init__(self, reference: Reference, t: Transcript, pad: int = 2000) -> None:
        self.reference = reference
        self.chrom = t.chrom
        self.start0 = max(0, min(s for s, _ in t.exons) - 1 - pad)
        self.end0 = min(reference.lengths[t.chrom], max(e for _, e in t.exons) + pad)
        self.seq = reference.fetch(t.chrom, self.start0, self.end0)

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        if chrom == self.chrom and self.start0 <= start0 and end0 <= self.end0:
            return self.seq[start0 - self.start0 : end0 - self.start0]
        return self.reference.fetch(chrom, start0, end0)


def nmd_eligible(t: Transcript, pos: int, cfg: Config) -> bool:
    exon = t.exon_of(pos)
    if exon < 0 or exon == len(t.exons) - 1:
        return False
    return t.tx_coord(pos) <= t.last_junction() - cfg.nmd_margin


def nonsense(t: Transcript, ref: GeneSequence, cfg: Config) -> Iterator[Candidate]:
    cds = t.cds_positions()
    for i in range(cfg.cds_start_exclude - cfg.cds_start_exclude % 3, len(cds) - 2, 3):
        codon_pos = cds[i : i + 3]
        if not nmd_eligible(t, codon_pos[0], cfg):
            continue
        genomic = [ref.fetch(t.chrom, p - 1, p) for p in codon_pos]
        codon = (
            "".join(genomic)
            if t.strand == "+"
            else "".join(b.translate(COMPLEMENT) for b in genomic)
        )
        if codon in STOPS:
            continue
        for j in range(3):
            for b in BASES:
                if b == codon[j]:
                    continue
                new = codon[:j] + b + codon[j + 1 :]
                if new in STOPS:
                    alt = b if t.strand == "+" else b.translate(COMPLEMENT)
                    yield Candidate(
                        "nonsense",
                        "-",
                        t.gene,
                        t.transcript_id,
                        t.chrom,
                        codon_pos[j],
                        genomic[j],
                        alt,
                        t.exon_of(codon_pos[j]) + 1,
                        len(t.exons),
                        f"codon {i // 3 + 1} {codon}>{new}",
                    )


def splice(t: Transcript, ref: GeneSequence, cfg: Config) -> Iterator[Candidate]:
    """Canonical donor +1/+2 and acceptor -1/-2 of junctions between NMD-eligible exons."""
    cds_set = set(t.cds_positions())
    for k in range(len(t.exons) - 1):
        up, down = t.exons[k], t.exons[k + 1]
        up_end = up[1] if t.strand == "+" else up[0]
        down_start = down[0] if t.strand == "+" else down[1]
        if up_end not in cds_set or down_start not in cds_set:
            continue  # junction outside the CDS
        if not nmd_eligible(t, up_end, cfg) or k + 1 == len(t.exons) - 1:
            continue
        step = 1 if t.strand == "+" else -1
        sites = [
            ("splice_donor", up_end + step),
            ("splice_donor", up_end + 2 * step),
            ("splice_acceptor", down_start - step),
            ("splice_acceptor", down_start - 2 * step),
        ]
        for kind, pos in sites:
            base = ref.fetch(t.chrom, pos - 1, pos)
            for b in BASES:
                if b != base:
                    yield Candidate(
                        kind,
                        "-",
                        t.gene,
                        t.transcript_id,
                        t.chrom,
                        pos,
                        base,
                        b,
                        k + 1 if kind == "splice_donor" else k + 2,
                        len(t.exons),
                        f"intron {k + 1}",
                    )


def frameshift(t: Transcript, ref: GeneSequence, cfg: Config) -> Iterator[Candidate]:
    """One deletion or insertion per CDS position and length bin; length and kind from the seed."""
    cds = t.cds_positions()
    for i in range(cfg.cds_start_exclude, len(cds)):
        p = cds[i]
        exon = t.exon_of(p)
        if exon < 0 or not nmd_eligible(t, p, cfg):
            continue
        for bin_name, all_lengths in cfg.frameshift_bins.items():
            lengths = [n for n in all_lengths if n % 3]
            pick = derive_int(cfg.seed, "synthetic_frameshift", f"{t.gene}|{p}|{bin_name}")
            kinds = cfg.frameshift_kinds
            length, kind = lengths[pick % len(lengths)], kinds[pick // 7 % len(kinds)]
            event = _frameshift_event(t, ref, cfg, p, exon, length, kind)
            if event is not None:
                anchor, r, a = event
                yield Candidate(
                    "frameshift",
                    bin_name,
                    t.gene,
                    t.transcript_id,
                    t.chrom,
                    anchor,
                    r,
                    a,
                    exon + 1,
                    len(t.exons),
                    f"{kind} {length} bp",
                )


def _frameshift_event(
    t: Transcript, ref: GeneSequence, cfg: Config, p: int, exon: int, length: int, kind: str
) -> tuple[int, str, str] | None:
    """VCF record (anchor, REF, ALT) of a deletion or insertion at CDS position ``p``, or None.

    The event must lie ``indel_exon_margin`` bases inside the exon and must not be shiftable
    (a repeat would let the same haplotype be written elsewhere).
    """
    s, e = t.exons[exon]
    m = cfg.indel_exon_margin
    if kind == "del":
        lo, hi = (p, p + length - 1) if t.strand == "+" else (p - length + 1, p)
        if lo - m < s or hi + m > e:
            return None
        anchor = lo - 1
        r = ref.fetch(t.chrom, anchor - 1, hi)
        a = r[0]
    else:
        anchor = p if t.strand == "+" else p - 1
        if anchor - m < s or anchor + 1 + m > e:
            return None
        key = f"{t.gene}|{p}|{length}|{kind}"
        ins = "".join(
            BASES[derive_int(cfg.seed, "synthetic_ins", f"{key}|{n}") % 4] for n in range(length)
        )
        r = ref.fetch(t.chrom, anchor - 1, anchor)
        a = r + ins
    if left_normalize(t.chrom, anchor, r, a, ref) != (anchor, r, a):
        return None
    if _shifts_right(ref, t.chrom, anchor, r, a):
        return None
    return anchor, r, a


def _shifts_right(ref: GeneSequence, chrom: str, pos: int, r: str, a: str) -> bool:
    """True when a pure indel could also be written one base to the right."""
    unit = r[1:] if len(r) > len(a) else a[1:]
    after = ref.fetch(chrom, pos - 1 + len(r), pos - 1 + len(r) + 1)
    return bool(unit) and after == unit[0]


def _group(kind: str) -> str:
    return "splice" if kind.startswith("splice") else kind


def ranked(cands: Iterable[Candidate], seed: int) -> list[Candidate]:
    return sorted(cands, key=lambda c: (derive_int(seed, "synthetic", c.key), c.key))


def generate(cfg: Config, out: Path) -> None:
    reference = Reference(cfg.reference)
    hi = haploinsufficient(cfg.clingen_dosage)
    genes = sorted(g for g in cfg.genes if g in hi and g not in cfg.exclude_genes)
    dropped_genes = sorted(set(cfg.genes) - set(genes))
    models = read_mane(cfg.mane_gtf, set(genes))
    missing = sorted(set(genes) - set(models))
    if missing:
        raise SystemExit(f"no MANE Select transcript for: {', '.join(missing)}")
    excl = Exclusions(cfg, reference)
    want = {k: int(n * cfg.oversample + 0.999) for k, n in cfg.targets.items()}
    pools: dict[str, list[Candidate]] = {k: [] for k in cfg.targets}
    for gene in genes:
        t = models[gene]
        seq = GeneSequence(reference, t)
        for gen in (nonsense, splice, frameshift):
            for c in gen(t, seq, cfg):
                group = (
                    c.bin
                    if c.kind == "frameshift"
                    else ("splice" if c.kind.startswith("splice") else c.kind)
                )
                if group in pools:
                    pools[group].append(c)
    chosen: list[Candidate] = []
    report = Counter[str]()
    for group, pool in pools.items():
        taken = 0
        for c in ranked(pool, cfg.seed):
            if taken >= want[group]:
                break
            why = excl.reason(c)
            if why is not None:
                report[f"{group}:{why}"] += 1
                continue
            chosen.append(c)
            taken += 1
        report[f"{group}:candidates"] = len(pool)
        report[f"{group}:kept"] = taken
    out.mkdir(parents=True, exist_ok=True)
    _write_vcf(out / "candidates.vcf", chosen, reference)
    _write_tsv(out / "candidates.tsv", chosen, cfg.seed)
    with (out / "generate_report.tsv").open("w", encoding="utf-8") as fh:
        fh.write(f"# genes\t{','.join(genes)}\n")
        fh.write(f"# genes dropped (not HI 3 or excluded)\t{','.join(dropped_genes)}\n")
        fh.write(f"# gnomad checked\t{cfg.gnomad_vcf is not None}\n")
        fh.write(f"# erepo checked\t{cfg.erepo_tsv is not None}\n")
        fh.write(f"# frameshift kinds\t{','.join(cfg.frameshift_kinds)}\n")
        fh.write(f"# deletion exclusion\t{cfg.deletion_exclusion}\n")
        for k in sorted(report):
            fh.write(f"{k}\t{report[k]}\n")
    reference.close()
    print(f"{len(chosen)} candidates -> {out}")


def _write_vcf(path: Path, cands: list[Candidate], reference: Reference) -> None:
    order = reference.order
    with path.open("w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for name, length in reference.contigs:
            fh.write(f"##contig=<ID={name},length={length}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c in sorted(cands, key=lambda c: (order[c.chrom], c.pos, c.ref, c.alt)):
            fh.write(f"{c.chrom}\t{c.pos}\t{c.key}\t{c.ref}\t{c.alt}\t.\tPASS\t.\n")


TSV_COLUMNS = [
    "key",
    "kind",
    "bin",
    "gene",
    "transcript",
    "chrom",
    "pos",
    "ref",
    "alt",
    "exon",
    "exons",
    "note",
    "rank",
]


def _write_tsv(path: Path, cands: list[Candidate], seed: int) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(TSV_COLUMNS)
        for c in ranked(cands, seed):
            w.writerow(
                [
                    c.key,
                    c.kind,
                    c.bin,
                    c.gene,
                    c.transcript,
                    c.chrom,
                    c.pos,
                    c.ref,
                    c.alt,
                    c.exon,
                    c.exons,
                    c.note,
                    derive_int(seed, "synthetic", c.key),
                ]
            )


# -- verify -------------------------------------------------------------------------------------


def read_vep(path: Path) -> dict[str, dict[str, dict[str, str]]]:
    """{variant key (ID column): {transcript id without version: CSQ fields}} from a VEP VCF."""
    fields: list[str] = []
    out: dict[str, dict[str, dict[str, str]]] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = line.split("Format: ")[1].split('"')[0].split("|")
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            per: dict[str, dict[str, str]] = {}
            for entry in info.get("CSQ", "").split(","):
                values = dict(zip(fields, entry.split("|"), strict=False))
                tx = values.get("Feature", "").split(".")[0]
                if tx:
                    per[tx] = values
            out[f[2]] = per
    return out


def verify(cfg: Config, out: Path, vep_path: Path) -> None:
    rows = list(csv.DictReader((out / "candidates.tsv").open(encoding="utf-8"), delimiter="\t"))
    vep = read_vep(vep_path)
    kept = Counter[str]()
    per_gene = Counter[str]()
    report = Counter[str]()
    placer = SamplePlacer(cfg)

    def group_of(r: dict[str, str]) -> str:
        return r["bin"] if r["kind"] == "frameshift" else _group(r["kind"])

    # round robin over the groups (rank order within each), so that no group uses up the room
    # (gene cap, spacing) of the others
    queues = {g: [r for r in rows if group_of(r) == g] for g in sorted(cfg.targets)}
    interleaved: list[dict[str, str]] = []
    while any(queues.values()):
        for q in queues.values():
            if q:
                interleaved.append(q.pop(0))
    for r in interleaved:
        group = group_of(r)
        if kept[group] >= cfg.targets[group]:
            continue
        csq = vep.get(r["key"], {}).get(r["transcript"].split(".")[0])
        if csq is None:
            report[f"{group}:no_vep_on_transcript"] += 1
            continue
        consequences = set(csq.get("Consequence", "").split("&"))
        if not consequences & EXPECTED_CONSEQUENCE[r["kind"]]:
            report[f"{group}:consequence_mismatch"] += 1
            continue
        if cfg.max_per_gene and per_gene[r["gene"]] >= cfg.max_per_gene:
            report[f"{group}:gene_cap"] += 1
            continue
        r = dict(
            r,
            hgvsc=csq.get("HGVSc", ""),
            hgvsp=csq.get("HGVSp", ""),
            consequence=csq.get("Consequence", ""),
        )
        if not placer.place(r):
            report[f"{group}:min_distance"] += 1
            continue
        kept[group] += 1
        per_gene[r["gene"]] += 1
    for group, n in cfg.targets.items():
        report[f"{group}:kept"] = kept[group]
        report[f"{group}:target"] = n
    by_sample = placer.out
    write_suite(cfg, by_sample)
    with (out / "synthetic_variants.tsv").open("w", newline="", encoding="utf-8") as fh:
        cols = [*TSV_COLUMNS, "hgvsc", "hgvsp", "consequence", "sample"]
        w = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        w.writeheader()
        for sample, rs in by_sample.items():
            for r in rs:
                w.writerow(dict(r, sample=sample))
    with (out / "verify_report.tsv").open("w", encoding="utf-8") as fh:
        for k in sorted(report):
            fh.write(f"{k}\t{report[k]}\n")
    short = [g for g in cfg.targets if report[f"{g}:kept"] < cfg.targets[g]]
    print(
        f"{sum(len(v) for v in by_sample.values())} variants -> {cfg.suite_out}"
        + (f"; short of target: {', '.join(short)}" if short else "")
    )


class SamplePlacer:
    """Keeps ``min_distance_bp`` between all variants (the framework's selection deduplication is
    suite-wide) and gives each variant to the sample with the fewest variants (rank order)."""

    def __init__(self, cfg: Config) -> None:
        self.min_distance = cfg.min_distance_bp
        self.spans: list[tuple[str, int, int]] = []
        self.out: dict[str, list[dict[str, str]]] = {s: [] for s in cfg.samples}

    def place(self, r: dict[str, str]) -> bool:
        s0, e0 = ref_span0(int(r["pos"]), r["ref"])
        if any(
            c == r["chrom"] and max(s0 - b, a - e0) < self.min_distance for c, a, b in self.spans
        ):
            return False
        self.spans.append((r["chrom"], s0, e0))
        sample = min(self.out, key=lambda s: (len(self.out[s]), s))
        self.out[sample].append(r)
        return True


def write_suite(cfg: Config, by_sample: dict[str, list[dict[str, str]]]) -> None:
    suite = yaml.safe_load(cfg.template.read_text(encoding="utf-8"))
    rules = suite.setdefault("selection", [])
    for sample, rows in by_sample.items():
        if not rows:
            continue
        variants = []
        for r in sorted(rows, key=lambda r: (r["chrom"], int(r["pos"]))):
            variants.append(
                {
                    "chrom": r["chrom"],
                    "pos": int(r["pos"]),
                    "ref": r["ref"],
                    "alt": r["alt"],
                    "gene": r["gene"],
                    "classification": "likely_pathogenic",
                    "acmg_criteria": ["PVS1", "PM2"],
                    "rationale": (
                        f"synthetic {r['kind']}{'' if r['bin'] == '-' else ' ' + r['bin'] + ' bp'} "
                        f"({r['note']}), {r['hgvsc']} {r['hgvsp']}, "
                        f"exon {r['exon']} of {r['exons']}, "
                        f"MANE Select {r['transcript']} (v{cfg.mane_release}); NMD predicted; "
                        "absent from ClinVar and eRepo"
                    ),
                }
            )
        rules.append(
            {
                "id": f"synthetic_{sample}",
                "category": "synthetic_novel",
                "source": "explicit",
                "variants": variants,
                "expected": {
                    "classification": ["pathogenic", "likely_pathogenic"],
                    "reportable": True,
                    "level": "both",
                    "basis": "synthetic LoF in an AD haploinsufficient SF gene: PVS1 + PM2",
                },
                "samples": [sample],
            }
        )
    header = (
        f"# Generated by scripts/generate_synthetic.py from {cfg.path.name} (seed {cfg.seed}); "
        "do not edit by hand.\n"
    )
    cfg.suite_out.write_text(
        header + yaml.safe_dump(suite, sort_keys=False, width=100), encoding="utf-8"
    )


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="step", required=True)
    for name in ("generate", "verify"):
        sp = sub.add_parser(name)
        sp.add_argument("--config", type=Path, required=True)
        sp.add_argument("--out", type=Path, required=True)
        if name == "generate":
            sp.add_argument("--gnomad", type=Path, help="gnomAD sites VCF (overrides the config)")
        if name == "verify":
            sp.add_argument("--vep", type=Path, required=True)
    a = ap.parse_args(argv)
    cfg = Config.load(a.config)
    if a.step == "generate":
        if a.gnomad is not None:
            cfg.gnomad_vcf = a.gnomad.resolve()
        generate(cfg, a.out)
    else:
        verify(cfg, a.out, a.vep)


if __name__ == "__main__":
    main()
