"""scripts/generate_synthetic.py on the fixture reference with hand-made gene models."""

from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pysam
import yaml

from tests.conftest import FIXTURES

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "generate_synthetic.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_synthetic", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["generate_synthetic"] = module
    spec.loader.exec_module(module)
    return module


gs = _load()

# GENEP on chr1 (+), GENEM on chr2 (-); exons and CDS 1-based closed
GTF = [
    ("chr1", "exon", 1001, 1200, "+", "GENEP", "TXP.1"),
    ("chr1", "exon", 2001, 2300, "+", "GENEP", "TXP.1"),
    ("chr1", "exon", 3001, 3400, "+", "GENEP", "TXP.1"),
    ("chr1", "CDS", 1101, 1200, "+", "GENEP", "TXP.1"),
    ("chr1", "CDS", 2001, 2300, "+", "GENEP", "TXP.1"),
    ("chr1", "CDS", 3001, 3150, "+", "GENEP", "TXP.1"),
    ("chr2", "exon", 7001, 7400, "-", "GENEM", "TXM.1"),
    ("chr2", "exon", 6001, 6300, "-", "GENEM", "TXM.1"),
    ("chr2", "exon", 5001, 5200, "-", "GENEM", "TXM.1"),
    ("chr2", "CDS", 7001, 7300, "-", "GENEM", "TXM.1"),
    ("chr2", "CDS", 6001, 6300, "-", "GENEM", "TXM.1"),
    ("chr2", "CDS", 5101, 5200, "-", "GENEM", "TXM.1"),
]


def _config(
    tmp_path: Path, clinvar: list[tuple[str, int, str, str]] | None = None, **extra: Any
) -> Any:
    gtf = tmp_path / "mane.gtf"
    attrs = 'gene_name "{}"; transcript_id "{}"; tag "MANE_Select";'
    gtf.write_text(
        "".join(
            f"{c}\tt\t{f}\t{s}\t{e}\t.\t{st}\t.\t{attrs.format(g, tx)}\n"
            for c, f, s, e, st, g, tx in GTF
        )
    )
    dosage = tmp_path / "dosage.tsv"
    dosage.write_text("#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHI\nGENEP\t1\t-\t-\t3\n"
                      "GENEM\t2\t-\t-\t3\nGENEX\t3\t-\t-\t1\n")  # fmt: skip
    header = pysam.VariantHeader()
    for name, length in (("1", 10000), ("2", 10000)):
        header.contigs.add(name, length=length)
    vcf = tmp_path / "clinvar.vcf.gz"
    with pysam.VariantFile(str(vcf), "wz", header=header) as out:
        for chrom, pos, ref, alt in clinvar or []:
            out.write(out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt)))
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    (tmp_path / "empty.bed").write_text("")
    template = tmp_path / "template.yaml"
    template.write_text(yaml.safe_dump({"schema_version": 1, "selection": []}))
    cfg = {
        "seed": 7, "reference": str(FIXTURES / "ref.fa"),
        "mane": {"gtf": "mane.gtf", "release": "test"}, "clingen_dosage": "dosage.tsv",
        "genes": ["GENEP", "GENEM", "GENEX"], "exclude_genes": [],
        "clinvar_vcf": "clinvar.vcf.gz", "avoid_strata": ["empty.bed"],
        "cds_start_exclude": 30, "nmd_margin": 50, "min_distance_bp": 100, "oversample": 1.0,
        "targets": {"nonsense": 5, "splice": 4, "1-5": 5, "6-10": 5},
        "frameshift_bins": {"1-5": [1, 2, 4, 5], "6-10": [7, 8, 10]},
        "samples": ["p1", "p2"], "template": "template.yaml", "suite_out": "suite.yaml",
    }  # fmt: skip
    cfg.update(extra)
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return gs.Config.load(path)


def _models(cfg: Any) -> dict[str, Any]:
    return gs.read_mane(cfg.mane_gtf, {"GENEP", "GENEM"})


def _seq(cfg: Any, t: Any) -> Any:
    from insilico_clinval.genome.reference import Reference

    return gs.GeneSequence(Reference(cfg.reference), t)


def _codon_after(t: Any, seq: Any, c: Any) -> str:
    """Codon (transcript orientation) containing the changed base, after the change."""
    cds = t.cds_positions()
    i = cds.index(c.pos)
    pos = cds[i - i % 3 : i - i % 3 + 3]
    bases = [c.alt if p == c.pos else seq.fetch(t.chrom, p - 1, p) for p in pos]
    codon = "".join(bases)
    return codon if t.strand == "+" else "".join(b.translate(gs.COMPLEMENT) for b in bases)


def test_gene_models_in_transcript_order(tmp_path: Path) -> None:
    models = _models(_config(tmp_path))
    m = models["GENEM"]
    assert m.exons[0] == (7001, 7400) and m.exons[-1] == (5001, 5200)
    assert m.cds_positions()[:2] == [7300, 7299]
    assert m.tx_coord(7400) == 1 and m.last_junction() == 400 + 300
    assert models["GENEP"].tx_coord(2001) == 201


def test_nonsense_on_both_strands(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    for t in _models(cfg).values():
        seq = _seq(cfg, t)
        cands = list(gs.nonsense(t, seq, cfg))
        assert cands
        for c in cands:
            assert _codon_after(t, seq, c) in gs.STOPS
            assert c.ref == seq.fetch(t.chrom, c.pos - 1, c.pos)
            assert t.exon_of(c.pos) != len(t.exons) - 1  # never in the last exon
            assert t.tx_coord(c.pos) <= t.last_junction() - cfg.nmd_margin
            assert t.cds_positions().index(c.pos) >= cfg.cds_start_exclude


def test_splice_sites_on_both_strands(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    models = _models(cfg)
    plus = {(c.kind, c.pos) for c in gs.splice(models["GENEP"], _seq(cfg, models["GENEP"]), cfg)}
    assert plus == {("splice_donor", 1201), ("splice_donor", 1202),
                    ("splice_acceptor", 2000), ("splice_acceptor", 1999)}  # fmt: skip
    minus = {(c.kind, c.pos) for c in gs.splice(models["GENEM"], _seq(cfg, models["GENEM"]), cfg)}
    assert minus == {("splice_donor", 7000), ("splice_donor", 6999),
                     ("splice_acceptor", 6301), ("splice_acceptor", 6302)}  # fmt: skip


def test_frameshifts(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    for t in _models(cfg).values():
        seq = _seq(cfg, t)
        cands = list(gs.frameshift(t, seq, cfg))
        assert {c.bin for c in cands} == {"1-5", "6-10"}
        for c in cands:
            assert abs(len(c.alt) - len(c.ref)) % 3 != 0
            assert c.ref[0] == c.alt[0]
            s0, e0 = c.span
            s, e = t.exons[c.exon - 1]
            assert s - 1 <= s0 and e0 <= e  # inside the exon
            assert gs.left_normalize(c.chrom, c.pos, c.ref, c.alt, seq) == (c.pos, c.ref, c.alt)
            assert c.ref == seq.fetch(c.chrom, c.pos - 1, c.pos - 1 + len(c.ref))


def test_generate_excludes_clinvar_and_is_deterministic(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    out_a, out_b = tmp_path / "a", tmp_path / "b"
    gs.generate(cfg, out_a)
    gs.generate(cfg, out_b)
    for name in ("candidates.vcf", "candidates.tsv", "generate_report.tsv"):
        assert (out_a / name).read_bytes() == (out_b / name).read_bytes()
    first = (out_a / "candidates.tsv").read_text().splitlines()[1].split("\t")
    chrom, pos, ref, alt = first[5], int(first[6]), first[7], first[8]
    blocked = tmp_path / "blocked"
    blocked.mkdir()
    cfg2 = _config(blocked, clinvar=[(chrom.removeprefix("chr"), pos, ref, alt)])
    gs.generate(cfg2, blocked / "out")
    keys = [
        line.split("\t")[0]
        for line in (blocked / "out" / "candidates.tsv").read_text().splitlines()
    ]
    assert first[0] not in keys
    assert "in_clinvar" in (blocked / "out" / "generate_report.tsv").read_text()
    assert "GENEX" in (out_a / "generate_report.tsv").read_text().splitlines()[1]  # HI 1: dropped


def test_verify_drops_mismatches_and_keeps_distance(tmp_path: Path) -> None:
    cfg = _config(tmp_path)
    out = tmp_path / "out"
    gs.generate(cfg, out)
    rows = [line.split("\t") for line in (out / "candidates.tsv").read_text().splitlines()[1:]]
    expected = {
        "nonsense": "stop_gained",
        "splice_donor": "splice_donor_variant",
        "splice_acceptor": "splice_acceptor_variant",
        "frameshift": "frameshift_variant",
    }
    wrong = rows[0][0]
    lines = ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
             'Ensembl VEP. Format: Allele|Consequence|Feature|HGVSc|HGVSp">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]  # fmt: skip
    for r in rows:
        csq = "missense_variant" if r[0] == wrong else expected[r[1]]
        lines.append(f"{r[5]}\t{r[6]}\t{r[0]}\t{r[7]}\t{r[8]}\t.\t.\tCSQ=A|{csq}|{r[4]}|c.1|p.1")
    vep = tmp_path / "vep.vcf"
    vep.write_text("\n".join(lines) + "\n")
    gs.verify(cfg, out, vep)
    report = (out / "verify_report.tsv").read_text()
    assert "consequence_mismatch\t1" in report
    suite = yaml.safe_load(cfg.suite_out.read_text())
    variants = {r["samples"][0]: r["variants"] for r in suite["selection"]}
    assert wrong not in {
        f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}" for vs in variants.values() for v in vs
    }
    spans = sorted(
        (v["chrom"], v["pos"], v["pos"] + len(v["ref"]) - 1) for vs in variants.values() for v in vs
    )  # suite-wide, like the framework
    for (c1, _, e1), (c2, s2, _) in itertools.pairwise(spans):
        assert c1 != c2 or s2 - 1 - e1 >= cfg.min_distance_bp
    sizes = [len(vs) for vs in variants.values()]
    assert max(sizes) - min(sizes) <= 1  # balanced
    assert all(v["acmg_criteria"] == ["PVS1", "PM2"] for vs in variants.values() for v in vs)


def test_deletions_only(tmp_path: Path) -> None:
    cfg = _config(tmp_path, frameshift_kinds=["del"])
    for t in _models(cfg).values():
        cands = list(gs.frameshift(t, _seq(cfg, t), cfg))
        assert cands and all(len(c.ref) > len(c.alt) for c in cands)


def test_deletion_exclusion_by_breakpoints(tmp_path: Path) -> None:
    """An interior ClinVar record excludes a deletion under ``span`` only; one at a breakpoint
    excludes it under both rules; insertions are unaffected by the rule."""
    from insilico_clinval.genome.reference import Reference

    pos, n = 2100, 10  # chr1, inside exon 2 of GENEP
    reference = Reference(FIXTURES / "ref.fa")
    ref = reference.fetch("chr1", pos - 1, pos + n)
    base = lambda p: reference.fetch("chr1", p - 1, p)  # noqa: E731
    deletion = gs.Candidate("frameshift", "6-10", "GENEP", "TXP.1", "chr1", pos, ref, ref[0],
                            2, 3, "del 10 bp")  # fmt: skip
    insertion = gs.Candidate("frameshift", "1-5", "GENEP", "TXP.1", "chr1", pos + 5,
                             base(pos + 5), base(pos + 5) + "A", 2, 3, "ins 1 bp")  # fmt: skip
    interior = ("1", pos + 5, base(pos + 5), "T" if base(pos + 5) != "T" else "G")
    edge = ("1", pos + n, base(pos + n), "T" if base(pos + n) != "T" else "G")  # last deleted base
    for records, rule, expected in [
        ([interior], "span", "in_clinvar"),
        ([interior], "breakpoints", None),
        ([edge], "breakpoints", "in_clinvar"),
    ]:
        d = tmp_path / f"{rule}_{records[0][1]}"
        d.mkdir()
        cfg = _config(d, clinvar=records, deletion_exclusion=rule)
        excl = gs.Exclusions(cfg, reference)
        assert excl.reason(deletion) == expected
        if records == [interior]:
            assert excl.reason(insertion) == "in_clinvar"  # anchor overlaps: as before
