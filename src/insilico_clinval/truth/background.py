"""Background genome (e.g. GIAB HG002 v4.2.1 GRCh38) merged into a sample's truth VCF.

The background VCF (bgzipped + tabix-indexed) is streamed contig by contig in reference order:

- records without an alternate allele in the background sample (``0/0``), with a missing
  genotype or with non-simple alleles are skipped and counted;
- records whose REF span overlaps a planted variant or lies closer than ``min_distance_bp``
  (``background.min_distance_bp``, default ``constraints.min_distance_bp``; same gap definition as
  selection deduplication) are removed and logged to
  ``<sample>.background_removed.tsv``;
- unphased heterozygous genotypes get a haplotype from
  ``derive_choice(seed, "bg_haplotype", f"{sample}|{chrom}:{pos}", 2)``; phased genotypes are kept
  and the background ``PS`` is dropped (relative phase of background blocks is not preserved);
- ALT alleles absent from the genotype are dropped; multi-allelic records are split later by
  ``bcftools norm -m -any``.

Optional filters (``BackgroundConfig``; all off by default, which keeps the v4.2.1 truth VCFs
unchanged): ``restrict_to_bed`` skips records whose REF span is not inside the BED,
``require_pass`` skips records with a FILTER other than PASS or missing, ``max_length`` skips
records with a genotyped allele longer than the limit (``event_length``), and
``spanning_deletion: ref`` reads ``*`` alleles as REF on their haplotype instead of skipping the
record (the deletion itself is carried by another record).

GIAB v4.2.1 covers chr1-22 only, so a male HG002 background has no chrX/chrY records; the
assembly-based v5.0q covers chrX/chrY too and needs the filters above (records outside its BED
can lack one haplotype).
"""

from __future__ import annotations

import bisect
import csv
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from insilico_clinval.genome.alleles import event_length, is_simple_allele, ref_span0
from insilico_clinval.genome.bed import IntervalSet
from insilico_clinval.genome.reference import Reference
from insilico_clinval.rng import derive_choice
from insilico_clinval.schemas.suite import BackgroundConfig
from insilico_clinval.schemas.variant import TestVariant
from insilico_clinval.truth import TruthError

REMOVED_COLUMNS = ["chrom", "pos", "ref", "alt", "gt", "near_variant_id", "gap_bp"]


@dataclass(frozen=True)
class BackgroundRecord:
    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    gt: tuple[int, ...]  # phased (diploid) or haploid


@dataclass
class BackgroundMerger:
    """Streams background records for one sample, removing those near planted variants."""

    config: BackgroundConfig
    vcf_path: Path
    bed_path: Path
    reference: Reference
    seed: int
    sample_id: str
    planted: list[TestVariant]
    min_distance_bp: int
    written: int = 0
    skipped: dict[str, int] = field(default_factory=dict)
    removed: list[list[str]] = field(default_factory=list)
    spanning_deletion_to_ref: int = 0

    def __post_init__(self) -> None:
        if not any(Path(f"{self.vcf_path}{ext}").is_file() for ext in (".tbi", ".csi")):
            raise TruthError(f"background VCF {self.vcf_path} must be bgzipped and tabix-indexed")
        self._vcf = pysam.VariantFile(str(self.vcf_path))
        samples = list(self._vcf.header.samples)
        self.column = self.config.sample or (samples[0] if samples else "")
        if self.column not in samples:
            raise TruthError(
                f"background VCF {self.vcf_path}: sample {self.column!r} not in {samples}"
            )
        self._contigs: dict[str, str] = {}
        for contig in self._vcf.header.contigs:
            name = str(contig)
            resolved = self.reference.resolve_contig(name)
            if resolved is not None:
                self._contigs[resolved] = name
        self._spans: dict[str, list[tuple[int, int, str]]] = {}
        for v in self.planted:
            start, end = ref_span0(v.pos, v.ref)
            self._spans.setdefault(v.chrom, []).append((start, end, v.variant_id))
        for spans in self._spans.values():
            spans.sort()
        self._starts = {c: [s for s, _, _ in spans] for c, spans in self._spans.items()}
        self._bed = IntervalSet(self.bed_path) if self.config.restrict_to_bed else None

    def _skip(self, reason: str) -> None:
        self.skipped[reason] = self.skipped.get(reason, 0) + 1

    def _near(self, chrom: str, start: int, end: int) -> tuple[str, int] | None:
        spans = self._spans.get(chrom)
        if not spans:
            return None
        # planted spans are short; scan the neighbourhood of the insertion point
        i = bisect.bisect_left(self._starts[chrom], start)
        for j in range(max(0, i - 2), min(len(spans), i + 2)):
            a_start, a_end, vid = spans[j]
            gap = max(start - a_end, a_start - end)
            if gap < self.min_distance_bp:
                return vid, gap
        return None

    def records(self, chrom: str) -> Iterator[BackgroundRecord]:
        name = self._contigs.get(chrom)
        if name is None:
            return
        cfg = self.config
        for rec in self._vcf.fetch(name):
            ref = str(rec.ref).upper()
            start, end = ref_span0(rec.pos, ref)
            if self._bed is not None and not self._bed.contains(chrom, start, end):
                self._skip("outside_bed")
                continue
            alts = tuple(a.upper() for a in (rec.alts or ()))
            spanning = {i for i, a in enumerate(alts, start=1) if a == "*"}
            checked = [a for i, a in enumerate(alts, start=1) if i not in spanning]
            if (
                not alts
                or (spanning and cfg.spanning_deletion == "skip")
                or not is_simple_allele(ref)
                or not all(map(is_simple_allele, checked))
            ):
                self._skip("non_simple_allele")
                continue
            if cfg.require_pass and any(f != "PASS" for f in rec.filter):
                self._skip("not_pass")
                continue
            call = rec.samples[self.column]
            gt = call.get("GT")
            if gt is None or any(a is None for a in gt):
                self._skip("missing_genotype")
                continue
            if spanning and any(a in spanning for a in gt):
                gt = tuple(0 if a in spanning else a for a in gt)
                self.spanning_deletion_to_ref += 1
            if all(a == 0 for a in gt):
                self._skip("no_alt_allele")
                continue
            if cfg.max_length is not None and any(
                event_length(ref, alts[a - 1]) > cfg.max_length for a in set(gt) if a > 0
            ):
                self._skip("too_long")
                continue
            near = self._near(chrom, start, end)
            if near is not None:
                sep = "|" if call.phased else "/"
                self.removed.append(
                    [
                        chrom,
                        str(rec.pos),
                        ref,
                        ",".join(alts),
                        sep.join(str(a) for a in gt),
                        near[0],
                        str(near[1]),
                    ]
                )
                continue
            alleles = tuple(int(a) for a in gt)
            if len(alleles) == 2 and not call.phased and alleles[0] != alleles[1]:
                key = f"{self.sample_id}|{chrom}:{rec.pos}"
                if derive_choice(self.seed, "bg_haplotype", key, 2) == 1:
                    alleles = (alleles[1], alleles[0])
                else:
                    alleles = (alleles[0], alleles[1])
            used = sorted({a for a in alleles if a > 0})
            remap = {0: 0, **{old: new for new, old in enumerate(used, start=1)}}
            self.written += 1
            yield BackgroundRecord(
                chrom=chrom,
                pos=rec.pos,
                ref=ref,
                alts=tuple(alts[a - 1] for a in used),
                gt=tuple(remap[a] for a in alleles),
            )

    def planted_outside_bed(self) -> list[str]:
        bed = IntervalSet(self.bed_path)
        return sorted(
            v.variant_id
            for v in self.planted
            if not bed.contains(v.chrom, *ref_span0(v.pos, v.ref))
        )

    def write_removed(self, path: Path) -> Path:
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(REMOVED_COLUMNS)
            writer.writerows(self.removed)
        return path

    def close(self) -> None:
        self._vcf.close()
