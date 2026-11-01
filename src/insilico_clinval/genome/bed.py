"""BED interval sets with overlap and containment queries.

Intervals are loaded into memory (merged, sorted), except for bgzipped BEDs with a tabix index,
which are queried on demand. BED contig names must match the reference.
"""

from __future__ import annotations

import bisect
import gzip
from pathlib import Path

import pysam


class IntervalSet:
    def __init__(self, path: Path):
        self.path = path
        self._tabix: pysam.TabixFile | None = None
        self._starts: dict[str, list[int]] = {}
        self._ends: dict[str, list[int]] = {}
        if path.suffix == ".gz" and Path(f"{path}.tbi").is_file():
            self._tabix = pysam.TabixFile(str(path))
            self._tabix_contigs = set(self._tabix.contigs)
        else:
            self._load()

    def _load(self) -> None:
        raw: dict[str, list[tuple[int, int]]] = {}
        opener = gzip.open if self.path.suffix == ".gz" else open
        with opener(self.path, "rt") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t")
                raw.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
        for chrom, intervals in raw.items():
            intervals.sort()
            starts: list[int] = []
            ends: list[int] = []
            for start, end in intervals:
                if starts and start <= ends[-1]:
                    ends[-1] = max(ends[-1], end)
                else:
                    starts.append(start)
                    ends.append(end)
            self._starts[chrom] = starts
            self._ends[chrom] = ends

    def _query(self, chrom: str, start0: int, end0: int) -> list[tuple[int, int]]:
        """Merged intervals overlapping the half-open span [start0, end0)."""
        if self._tabix is not None:
            if chrom not in self._tabix_contigs:
                return []
            hits = sorted(
                (int(f[1]), int(f[2]))
                for f in (row.split("\t") for row in self._tabix.fetch(chrom, start0, end0))
            )
            merged: list[tuple[int, int]] = []
            for s, e in hits:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            return merged
        starts = self._starts.get(chrom)
        if not starts:
            return []
        ends = self._ends[chrom]
        i = bisect.bisect_right(starts, start0) - 1
        i = max(i, 0)
        out: list[tuple[int, int]] = []
        while i < len(starts) and starts[i] < end0:
            if ends[i] > start0:
                out.append((starts[i], ends[i]))
            i += 1
        return out

    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        return bool(self._query(chrom, start0, end0))

    def contains(self, chrom: str, start0: int, end0: int) -> bool:
        return any(s <= start0 and end0 <= e for s, e in self._query(chrom, start0, end0))
