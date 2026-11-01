"""Simulated reference segments and their ploidy.

Region mode: every planted variant ± ``flank_bp`` (clipped to the contig), plus optional extra
regions (e.g. paralogous loci), merged per contig. Genome mode: whole contigs (default chr1-22,
chrX, chrY). In both modes, for male samples chrX and chrY are split at the pseudo-autosomal
regions: chrX PARs are diploid (they stand for the chrX and chrY copies), the rest of chrX and
chrY is haploid, and chrY PARs are left out. Female samples have no chrY segments.
"""

from __future__ import annotations

import bisect
from collections.abc import Iterable
from pathlib import Path

from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.enums import Sex
from insilico_clinval.schemas.simulation import SegmentRecord
from insilico_clinval.schemas.suite import Region

DEFAULT_GENOME_CONTIGS: tuple[str, ...] = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY")

Interval = tuple[str, int, int]  # chrom, 1-based start, 1-based end (closed)


def merge_intervals(intervals: Iterable[Interval], reference: Reference) -> list[Interval]:
    """Merge overlapping or adjacent intervals; sorted in FASTA-index order."""
    ordered = sorted(intervals, key=lambda iv: (reference.order[iv[0]], iv[1], iv[2]))
    merged: list[Interval] = []
    for chrom, start, end in ordered:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def read_bed_intervals(path: Path, reference: Reference) -> list[Interval]:
    """BED (0-based half-open) -> 1-based closed intervals; contigs must exist in the reference."""
    out: list[Interval] = []
    with path.open() as fh:
        for n, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = reference.resolve_contig(fields[0])
            if chrom is None:
                raise ValueError(f"{path}:{n}: contig {fields[0]} not in the reference")
            start0, end0 = int(fields[1]), int(fields[2])
            if end0 <= start0:
                raise ValueError(f"{path}:{n}: empty interval")
            out.append((chrom, start0 + 1, min(end0, reference.lengths[chrom])))
    return out


def _split_by_ploidy(
    interval: Interval, sex: Sex, pars: tuple[Region, ...]
) -> list[tuple[Interval, int]]:
    """Split one interval at PAR boundaries; returns (piece, ploidy); chrY PARs are dropped."""
    chrom, start, end = interval
    if sex == Sex.FEMALE:
        return [] if chrom == "chrY" else [(interval, 2)]
    if chrom not in ("chrX", "chrY"):
        return [(interval, 2)]
    cuts = sorted(
        (max(r.start, start), min(r.end, end))
        for r in pars
        if r.chrom == chrom and r.start <= end and r.end >= start
    )
    pieces: list[tuple[Interval, int]] = []
    cursor = start
    for par_start, par_end in cuts:
        if cursor < par_start:
            pieces.append(((chrom, cursor, par_start - 1), 1))
        if chrom == "chrX":
            pieces.append(((chrom, par_start, par_end), 2))
        cursor = par_end + 1
    if cursor <= end:
        pieces.append(((chrom, cursor, end), 1))
    return pieces


def segment_name(chrom: str, start: int, end: int) -> str:
    return f"{chrom}_{start}_{end}"


def build_segments(
    intervals: Iterable[Interval], reference: Reference, sex: Sex, pars: tuple[Region, ...]
) -> list[SegmentRecord]:
    segments: list[SegmentRecord] = []
    for interval in merge_intervals(intervals, reference):
        for (chrom, start, end), ploidy in _split_by_ploidy(interval, sex, pars):
            segments.append(
                SegmentRecord(
                    name=segment_name(chrom, start, end),
                    chrom=chrom,
                    start=start,
                    end=end,
                    ploidy=1 if ploidy == 1 else 2,
                )
            )
    return segments


def region_intervals(
    variants: Iterable[tuple[str, int, str]],
    reference: Reference,
    flank_bp: int,
    extra: Iterable[Interval] = (),
) -> list[Interval]:
    """Planted variants (chrom, pos, ref) ± flank, clipped to the contig, plus extra intervals."""
    out: list[Interval] = list(extra)
    for chrom, pos, ref in variants:
        length = reference.lengths[chrom]
        out.append((chrom, max(1, pos - flank_bp), min(length, pos + len(ref) - 1 + flank_bp)))
    return out


def genome_intervals(reference: Reference, contigs: list[str] | None) -> list[Interval]:
    if contigs is None:
        names = [c for c in DEFAULT_GENOME_CONTIGS if c in reference.lengths]
    else:
        missing = [c for c in contigs if c not in reference.lengths]
        if missing:
            raise ValueError(f"simulation.genome_contigs not in the reference: {missing}")
        names = list(contigs)
    return [(c, 1, reference.lengths[c]) for c in names]


class SegmentIndex:
    """Lookup of the segment overlapping / containing a genome interval."""

    def __init__(self, segments: list[SegmentRecord]):
        self._by_chrom: dict[str, list[SegmentRecord]] = {}
        for seg in sorted(segments, key=lambda s: (s.chrom, s.start)):
            self._by_chrom.setdefault(seg.chrom, []).append(seg)
        self._starts = {c: [s.start for s in segs] for c, segs in self._by_chrom.items()}

    def overlapping(self, chrom: str, start: int, end: int) -> SegmentRecord | None:
        """The segment overlapping the closed interval [start, end] (segments are disjoint)."""
        segs = self._by_chrom.get(chrom)
        if not segs:
            return None
        i = bisect.bisect_right(self._starts[chrom], end) - 1
        if i >= 0 and segs[i].end >= start:
            return segs[i]
        return None

    def contains(self, chrom: str, start: int, end: int) -> bool:
        seg = self.overlapping(chrom, start, end)
        return seg is not None and seg.start <= start and end <= seg.end
