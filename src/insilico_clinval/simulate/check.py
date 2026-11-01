"""Planting-completeness check from Mason's true alignments.

For each planted variant a reference window is defined: the REF span, extended to the right by
the distance an indel can be shifted inside a repeat (Mason's realignment may place the gap
anywhere in it) and padded on both sides by ``pad + max(len(REF), len(ALT))`` bases. The length
term matters: Mason realigns each read to the reference with edit-distance scoring and free end
gaps, so an event closer to a read end than its own length is represented by a shifted start
instead of a gap; such reads do not reach past the padded window and are ignored. Every primary
read that spans the whole window is classified first by its sequence: a read that contains the
whole ALT window sequence but not the REF one (or the reverse) is ``alt`` (``ref``). This does not
depend on Mason's CIGAR, which for a tandem duplication inside a repeat can scatter the inserted
bases over many small gaps. Otherwise the read bases aligned between the window ends are compared
with the REF and ALT window sequences (exact match first, then edit distance; ties count as
``other``). Counts are kept per Mason haplotype (``oH`` tag).

Indels longer than ``LONG_INDEL_BP`` (outside the reference implementation's 1-20 bp) have a
window longer than a read, so no read spans it; indels that can shift in a repeat have a wider
window, and Mason's shifted starts leave few reads spanning it. For both, a read that overlaps the
window without spanning it is classified by k-mers of the local haplotypes instead. The local
haplotype of an allele
is ``JUNCTION_BP`` bases of reference before it, the allele, the reference over which the indel
can shift in a repeat, and ``JUNCTION_BP`` bases after that; its marker k-mers (length
``2 * JUNCTION_BP``) are those absent from the other allele's local haplotype, so they work for
deletions, insertions and tandem duplications alike. A read that contains an ALT marker and no
REF marker is ``alt`` (and vice versa); other reads are not counted.

A variant is ``ok`` when every carrying haplotype has at least ``min_alt_reads`` ALT reads,
``missing_on_haplotype`` when a carrying haplotype has spanning reads but too few ALT reads, and
``no_reads`` when no read spans the window on some carrying haplotype.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from pathlib import Path

import pysam

from insilico_clinval.genome.alleles import event_length
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import Haplotype, PlantingStatus, Zygosity
from insilico_clinval.schemas.sample import SampleVariant
from insilico_clinval.schemas.simulation import HaplotypeReads, PlantingCheck, VariantPlanting
from insilico_clinval.simprep.prepare import TruthRecord

WINDOW_PAD = 3
MAX_SHIFT = 1000
LONG_INDEL_BP = 20
JUNCTION_BP = 20


def right_shift(reference: Reference, chrom: str, pos: int, ref: str, alt: str) -> int:
    """How far a pure insertion/deletion can move to the right without changing the haplotype."""
    if len(ref) == len(alt) or (len(ref) > 1 and len(alt) > 1) or ref[0] != alt[0]:
        return 0
    unit = ref[1:] if len(ref) > len(alt) else alt[1:]
    start0 = pos  # first base after the anchor (0-based)
    following = reference.fetch(
        chrom, start0, min(start0 + len(unit) + MAX_SHIFT, reference.lengths[chrom])
    )
    shift = 0
    if len(ref) > len(alt):
        # deletion of ref[1:]: shiftable while the base after the deleted run equals its first base
        k = len(unit)
        while k + shift < len(following) and following[shift] == following[k + shift]:
            shift += 1
    else:
        while shift < len(following) and following[shift] == unit[shift % len(unit)]:
            shift += 1
    return shift


def edit_distance(a: str, b: str) -> int:
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        current = [i]
        for j, cb in enumerate(b, start=1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (ca != cb)))
        previous = current
    return previous[-1]


def classify_by_sequence(sequence: str | None, ref_window: str, alt_window: str) -> str | None:
    """``alt``/``ref`` when the read contains exactly one of the two window sequences."""
    if sequence is None:
        return None
    has_alt, has_ref = alt_window in sequence, ref_window in sequence
    if has_alt == has_ref:
        return None
    return "alt" if has_alt else "ref"


def marker_kmers(own: str, other: str, k: int = 2 * JUNCTION_BP) -> set[str]:
    """k-mers of the local haplotype ``own`` that do not occur in the local haplotype ``other``."""
    return {own[i : i + k] for i in range(len(own) - k + 1)} - {
        other[i : i + k] for i in range(len(other) - k + 1)
    }


def classify_by_markers(
    sequence: str | None, ref_points: set[str], alt_points: set[str]
) -> str | None:
    """``alt``/``ref`` when the read contains marker k-mers of exactly one allele."""
    if sequence is None:
        return None
    has_alt = any(b in sequence for b in alt_points)
    has_ref = any(b in sequence for b in ref_points)
    if has_alt == has_ref:
        return None
    return "alt" if has_alt else "ref"


def classify_read_allele(observed: str, ref_window: str, alt_window: str) -> str:
    if observed == alt_window:
        return "alt"
    if observed == ref_window:
        return "ref"
    d_alt, d_ref = edit_distance(observed, alt_window), edit_distance(observed, ref_window)
    if d_alt < d_ref:
        return "alt"
    if d_ref < d_alt:
        return "ref"
    return "other"


def window_bases(read: pysam.AlignedSegment, start0: int, end0: int) -> str | None:
    """Read bases aligned to the reference window [start0, end0), insertions included."""
    q_start = q_end = None
    for q, r in read.get_aligned_pairs():
        if r == start0 and q is not None:
            q_start = q
        if r == end0 - 1 and q is not None:
            q_end = q
    if q_start is None or q_end is None or q_end < q_start:
        return None
    sequence = read.query_sequence
    return None if sequence is None else sequence[q_start : q_end + 1]


def carrying_haplotypes(sv: SampleVariant) -> list[int]:
    if sv.zygosity in (Zygosity.HEMI,) or sv.haplotype == Haplotype.ONE:
        return [1]
    if sv.haplotype == Haplotype.TWO:
        return [2]
    return [1, 2]


def expected_vaf(sv: SampleVariant) -> float:
    return 0.5 if sv.zygosity == Zygosity.HET else 1.0


def variant_window(
    reference: Reference, chrom: str, pos: int, ref: str, alt: str, pad: int = WINDOW_PAD
) -> tuple[int, int]:
    """Half-open 0-based reference window of a variant (see module docstring)."""
    shift = right_shift(reference, chrom, pos, ref, alt)
    margin = pad + max(len(ref), len(alt))
    start0 = max(0, pos - 1 - margin)
    end0 = min(reference.lengths[chrom], pos - 1 + len(ref) + shift + margin)
    return start0, end0


def check_variant(
    bam: pysam.AlignmentFile,
    reference: Reference,
    record: TruthRecord,
    sv: SampleVariant,
    *,
    min_alt_reads: int,
    pad: int = WINDOW_PAD,
) -> VariantPlanting:
    chrom, pos, ref, alt = record.chrom, record.pos, record.ref, record.alt
    start0, end0 = variant_window(reference, chrom, pos, ref, alt, pad)
    ref_window = reference.fetch(chrom, start0, end0)
    left = ref_window[: pos - 1 - start0]
    right = ref_window[pos - 1 - start0 + len(ref) :]
    alt_window = left + alt + right
    shift = right_shift(reference, chrom, pos, ref, alt)
    by_markers = len(ref) != len(alt) and (event_length(ref, alt) > LONG_INDEL_BP or shift > 0)
    if by_markers:
        after0 = pos - 1 + len(ref)
        middle = reference.fetch(chrom, after0, after0 + shift)
        flank = reference.fetch(chrom, after0 + shift, after0 + shift + JUNCTION_BP)
        lk = reference.fetch(chrom, max(0, pos - 1 - JUNCTION_BP), pos - 1)
        ref_local, alt_local = lk + ref + middle + flank, lk + alt + middle + flank
        ref_points = marker_kmers(ref_local, alt_local)
        alt_points = marker_kmers(alt_local, ref_local)
    counts: dict[int, Counter[str]] = {}
    for read in bam.fetch(chrom, start0, end0):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start > start0 or (read.reference_end or 0) < end0:
            if not by_markers:
                continue
            by_points = classify_by_markers(read.query_sequence, ref_points, alt_points)
            if by_points is not None:
                hap = int(read.get_tag("oH")) if read.has_tag("oH") else 0
                counts.setdefault(hap, Counter())[by_points] += 1
            continue
        allele = classify_by_sequence(read.query_sequence, ref_window, alt_window)
        if allele is None:
            observed = window_bases(read, start0, end0)
            allele = (
                "other"
                if observed is None
                else classify_read_allele(observed, ref_window, alt_window)
            )
        hap = int(read.get_tag("oH")) if read.has_tag("oH") else 0
        counts.setdefault(hap, Counter())[allele] += 1
    carrying = carrying_haplotypes(sv)
    by_hap = [
        HaplotypeReads(
            haplotype=h,
            reads=sum(c.values()),
            alt=c["alt"],
            ref=c["ref"],
            other=c["other"],
        )
        for h, c in sorted(counts.items())
    ]
    total: Counter[str] = Counter()
    for c in counts.values():
        total.update(c)
    spanning = sum(total.values())
    status = PlantingStatus.OK
    for h in carrying:
        c = counts.get(h, Counter())
        if sum(c.values()) == 0:
            status = PlantingStatus.NO_READS
            break
        if c["alt"] < min_alt_reads:
            status = PlantingStatus.MISSING_ON_HAPLOTYPE
    return VariantPlanting(
        variant_id=sv.variant_id,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        zygosity=sv.zygosity,
        carrying_haplotypes=carrying,
        expected_vaf=expected_vaf(sv),
        spanning_reads=spanning,
        alt_reads=total["alt"],
        ref_reads=total["ref"],
        other_reads=total["other"],
        vaf=round(total["alt"] / spanning, 4) if spanning else None,
        by_haplotype=by_hap,
        status=status,
    )


def planting_check(
    bam_path: Path,
    reference: Reference,
    sample_id: str,
    planted: Iterable[TruthRecord],
    sample_variants: Iterable[SampleVariant],
    *,
    min_alt_reads: int,
) -> PlantingCheck:
    by_id = {sv.variant_id: sv for sv in sample_variants}
    results: list[VariantPlanting] = []
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for record in planted:
            assert record.variant_id is not None
            results.append(
                check_variant(
                    bam, reference, record, by_id[record.variant_id], min_alt_reads=min_alt_reads
                )
            )
    results.sort(key=lambda v: v.variant_id)
    return PlantingCheck(
        schema_version=SCHEMA_VERSION,
        sample_id=sample_id,
        min_alt_reads=min_alt_reads,
        window_pad=WINDOW_PAD,
        complete=all(v.status == PlantingStatus.OK for v in results),
        variants=results,
    )


def write_check_tsv(check: PlantingCheck, path: Path) -> Path:
    columns = [
        "variant_id", "chrom", "pos", "ref", "alt", "zygosity", "expected_vaf", "vaf",
        "spanning_reads", "alt_reads", "ref_reads", "other_reads", "alt_by_haplotype", "status",
    ]  # fmt: skip
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for v in check.variants:
            by_hap = ",".join(f"{h.haplotype}:{h.alt}/{h.reads}" for h in v.by_haplotype)
            row = [
                v.variant_id, v.chrom, str(v.pos), v.ref, v.alt, v.zygosity.value,
                str(v.expected_vaf), "" if v.vaf is None else str(v.vaf), str(v.spanning_reads),
                str(v.alt_reads), str(v.ref_reads), str(v.other_reads), by_hap, v.status.value,
            ]  # fmt: skip
            fh.write("\t".join(row) + "\n")
    return path
