"""Evidence for one variant from the intermediate outputs of the system under test.

- reads carrying the allele in the pipeline's alignment (same window and allele test as the
  planting check, so the counts are comparable with Mason's true alignments);
- records overlapping the variant in the raw VCF (PASS or not);
- the hap.py decision for a truth or query allele;
- the finding of the tertiary tool for an allele.

Alleles are compared after left-normalization against the reference, so representation
differences between the truth, the caller and the tertiary tool do not break the matching.
"""

from __future__ import annotations

import bisect
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from insilico_clinval.genome.alleles import left_normalize
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.finding import CanonicalFinding
from insilico_clinval.simulate.check import classify_read_allele, variant_window, window_bases

Allele = tuple[str, int, str, str]  # chrom, pos (1-based), ref, alt


def normalized(reference: Reference, allele: Allele) -> Allele:
    chrom, pos, ref, alt = allele
    if ref.upper() == alt.upper() or not ref or not alt or "*" in alt or "<" in alt:
        return chrom, pos, ref.upper(), alt.upper()
    pos, ref, alt = left_normalize(chrom, pos, ref, alt, reference)
    return chrom, pos, ref, alt


# -- alignment ---------------------------------------------------------------------------------


@dataclass(frozen=True)
class AlignmentCounts:
    spanning: int
    alt: int
    alt_mapq0: int


def alignment_counts(
    aln: pysam.AlignmentFile, reference: Reference, allele: Allele
) -> AlignmentCounts:
    """Primary reads spanning the variant window, and those carrying the ALT allele."""
    chrom, pos, ref, alt = allele
    start0, end0 = variant_window(reference, chrom, pos, ref, alt)
    ref_window = reference.fetch(chrom, start0, end0)
    alt_window = ref_window[: pos - 1 - start0] + alt + ref_window[pos - 1 - start0 + len(ref) :]
    spanning = alt_n = alt_q0 = 0
    for read in aln.fetch(chrom, start0, end0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        if read.reference_start > start0 or (read.reference_end or 0) < end0:
            continue
        spanning += 1
        observed = window_bases(read, start0, end0)
        if observed is not None and classify_read_allele(observed, ref_window, alt_window) == "alt":
            alt_n += 1
            alt_q0 += read.mapping_quality == 0
    return AlignmentCounts(spanning=spanning, alt=alt_n, alt_mapq0=alt_q0)


# -- raw VCF -----------------------------------------------------------------------------------


@dataclass(frozen=True)
class CallRecord:
    allele: Allele  # normalized
    passed: bool
    filters: str
    genotype: str


def calls_near(
    vcf: pysam.VariantFile, reference: Reference, chrom: str, start0: int, end0: int
) -> list[CallRecord]:
    """Records of the first sample overlapping [start0, end0), one per ALT allele."""
    out: list[CallRecord] = []
    if chrom not in vcf.header.contigs:
        return out
    for rec in vcf.fetch(chrom, start0, end0):
        filters = ";".join(str(k) for k in rec.filter) or "."
        passed = filters in ("PASS", ".")
        gt = rec.samples[0].get("GT") if len(rec.samples) else None
        genotype = "/".join("." if a is None else str(a) for a in gt) if gt else "."
        for alt in rec.alts or ():
            allele = normalized(reference, (rec.chrom, rec.pos, rec.ref or "", alt))
            out.append(CallRecord(allele=allele, passed=passed, filters=filters, genotype=genotype))
    return out


# -- hap.py annotated VCF ----------------------------------------------------------------------


@dataclass(frozen=True)
class HappyCall:
    allele: Allele  # normalized
    decision: str  # BD: TP, FP, FN, N, UNK
    kind: str  # BK: gm, am, lm, .


def _sample_value(rec: pysam.VariantRecord, sample: str, key: str) -> str:
    try:
        value = rec.samples[sample].get(key)
    except KeyError:
        return "."
    if value is None:
        return "."
    return str(value[0] if isinstance(value, tuple) else value)


def happy_calls(
    vcf: pysam.VariantFile,
    reference: Reference,
    sample: str,
    chrom: str,
    start0: int,
    end0: int,
) -> list[HappyCall]:
    """Decisions of ``sample`` (TRUTH or QUERY) for alleles overlapping [start0, end0)."""
    out: list[HappyCall] = []
    if chrom not in vcf.header.contigs:
        return out
    for rec in vcf.fetch(chrom, start0, end0):
        decision = _sample_value(rec, sample, "BD")
        if decision in (".", ""):
            continue
        kind = _sample_value(rec, sample, "BK")
        alleles = rec.alleles or ()
        gt = rec.samples[sample].get("GT") or ()
        alt_indexes = sorted({a for a in gt if a not in (None, 0)}) or range(1, len(alleles))
        for i in alt_indexes:
            allele = normalized(reference, (rec.chrom, rec.pos, rec.ref or "", alleles[i]))
            out.append(HappyCall(allele=allele, decision=decision, kind=kind))
    return out


def pick(calls: Sequence[HappyCall], allele: Allele) -> tuple[HappyCall | None, str]:
    """The call of exactly this allele, else (as a fallback) the first call in the window."""
    for c in calls:
        if c.allele == allele:
            return c, "exact"
    return (calls[0], "window") if calls else (None, "none")


# -- findings ----------------------------------------------------------------------------------


class FindingIndex:
    """Findings by position; matching normalizes only the findings near the queried allele."""

    def __init__(self, findings: Iterable[CanonicalFinding], reference: Reference):
        self.reference = reference
        self._by_chrom: dict[str, list[CanonicalFinding]] = {}
        for f in findings:
            self._by_chrom.setdefault(f.chrom, []).append(f)
        for items in self._by_chrom.values():
            items.sort(key=lambda f: f.pos)
        self._pos = {c: [f.pos for f in items] for c, items in self._by_chrom.items()}

    def find(self, allele: Allele, slack: int = 100) -> CanonicalFinding | None:
        chrom, pos, ref, _ = allele
        positions = self._pos.get(chrom, [])
        lo = bisect.bisect_left(positions, pos - slack)
        hi = bisect.bisect_right(positions, pos + len(ref) + slack)
        for f in self._by_chrom[chrom][lo:hi] if positions else []:
            if normalized(self.reference, (f.chrom, f.pos, f.ref, f.alt)) == allele:
                return f
        return None


def open_alignment(path: Path, reference_fasta: Path) -> pysam.AlignmentFile:
    return pysam.AlignmentFile(str(path), reference_filename=str(reference_fasta))
