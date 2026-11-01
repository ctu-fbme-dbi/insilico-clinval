"""Regions in which calls are scored against the truth.

The simulated genome is the reference plus the truth VCF, so the truth is complete wherever reads
were simulated: calls are scored in the sample's segments, trimmed at both ends so that the
coverage drop near a segment edge does not count as missed calls. This includes regions outside the
background benchmark's confident BED and chromosomes without background (chrX of HG002 v4.2.1),
where every call except a planted variant is a false positive. The confident BED is reported as a
stratum instead, for comparison with real data. A window around every planted variant is always
added, so a planted variant near a segment end is still scored.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from pathlib import Path

from insilico_clinval.genome.reference import Reference

Interval = tuple[str, int, int]  # chrom, 0-based start, end (half-open)
# hap.py 0.3.15 `quantify` crashes (segmentation fault) on a regions BED with whole-chromosome
# intervals; pieces of at most 1 Mb load fine and score the same bases.
MAX_INTERVAL_BP = 1_000_000


def read_segments(path: Path) -> list[Interval]:
    """Simulated segments of a sample (`<sample>.sim/segments.tsv`, 1-based closed)."""
    with path.open(newline="") as fh:
        return [
            (r["chrom"], int(r["start"]) - 1, int(r["end"]))
            for r in csv.DictReader(fh, delimiter="\t")
        ]


def planted_window(chrom: str, pos: int, ref: str, alt: str, pad: int) -> Interval:
    margin = pad + max(len(ref), len(alt))
    return chrom, max(0, pos - 1 - margin), pos - 1 + len(ref) + margin


def merge(intervals: Iterable[Interval], reference: Reference) -> list[Interval]:
    """Merged intervals sorted by reference contig order (contigs absent from it are dropped)."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e in intervals:
        if chrom in reference.order and e > s:
            by_chrom.setdefault(chrom, []).append((s, min(e, reference.lengths[chrom])))
    out: list[Interval] = []
    for chrom in sorted(by_chrom, key=lambda c: reference.order[c]):
        spans = sorted(by_chrom[chrom])
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append((chrom, cur_s, cur_e))
    return out


def evaluation_regions(
    segments: Sequence[Interval],
    windows: Sequence[Interval],
    reference: Reference,
    *,
    edge_bp: int,
) -> list[Interval]:
    """Trimmed segments plus planted windows, merged and sorted."""
    trimmed = ((c, s + edge_bp, e - edge_bp) for c, s, e in segments)
    return merge([*trimmed, *windows], reference)


def split(intervals: Iterable[Interval], max_bp: int = MAX_INTERVAL_BP) -> list[Interval]:
    """Cut every interval into consecutive pieces of at most ``max_bp``."""
    return [
        (chrom, s, min(s + max_bp, e))
        for chrom, start, e in intervals
        for s in range(start, e, max_bp)
    ]


def write_bed(intervals: Iterable[Interval], path: Path) -> Path:
    path.write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in intervals), encoding="utf-8")
    return path
