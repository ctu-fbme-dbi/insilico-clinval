"""Merging Mason outputs of the ploidy groups: FASTQ (gzip) and true alignments (genome coords).

FASTQ files are concatenated group by group into one deterministic gzip stream, compressed with
ISA-L (``isal.igzip_threaded``, the library behind ``igzip``). The header carries no file name and
``mtime=0``; the stream is cut into fixed blocks compressed independently, so the bytes do not
depend on the thread count (they do depend on the ISA-L level and the ``isal`` version). Mason's
alignments refer to segment contigs (``chr1_1001_4000``); they are lifted back to genome
coordinates (contig, position, mate contig and position) with the segment map and written as an
unsorted BAM whose header lists every contig of the reference FASTA. BGZF compression of that BAM
is multi-threaded; its output does not depend on the thread count either.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import IO

import pysam
from isal import igzip_threaded

from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.simulation import SegmentRecord

# ISA-L levels are 0-3; 3 is still ~10x faster than zlib level 6 and about 5 % larger.
ISAL_LEVEL = 3
# Also the gzip block size: copying in whole blocks keeps the block boundaries fixed.
_COPY_CHUNK = 1 << 20
LIFT_LOG_EVERY = 10_000_000

log = logging.getLogger(__name__)


class FastqGzConcat:
    """A deterministic gzip file that plain FASTQ files are appended to, one after another."""

    def __init__(self, dst: Path, *, threads: int = 1):
        self.path = dst
        self._gz: IO[bytes] = igzip_threaded.open(
            dst, "wb", compresslevel=ISAL_LEVEL, threads=threads, block_size=_COPY_CHUNK
        )

    def append(self, src: Path) -> None:
        with src.open("rb") as fh:
            shutil.copyfileobj(fh, self._gz, _COPY_CHUNK)

    def close(self) -> None:
        self._gz.close()

    def __enter__(self) -> FastqGzConcat:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def merge_fastq_gz(sources: Sequence[Path], dst: Path, *, threads: int = 1) -> Path:
    """Concatenate plain FASTQ files into one deterministic gzip file."""
    with FastqGzConcat(dst, threads=threads) as gz:
        for src in sources:
            gz.append(src)
    return dst


def genome_header(reference: Reference) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": name, "LN": length} for name, length in reference.contigs],
        }
    )


class AlignmentLifter:
    """Unsorted BAM in genome coordinates that Mason's alignment files are lifted into.

    ``threads`` is the number of BGZF (de)compression threads. With ``progress`` (a label, e.g.
    the sample id) the running count is logged every ``log_every`` alignments, as a percentage of
    ``expected`` when given.
    """

    def __init__(
        self,
        segments: Sequence[SegmentRecord],
        reference: Reference,
        dst: Path,
        *,
        threads: int = 1,
        expected: int | None = None,
        progress: str | None = None,
        log_every: int = LIFT_LOG_EVERY,
    ):
        self._by_name = {s.name: s for s in segments}
        self._header = genome_header(reference)
        self._threads = threads
        self._expected = expected
        self._progress = progress
        self._log_every = log_every
        self.count = 0
        self._out = pysam.AlignmentFile(str(dst), "wb", header=self._header, threads=threads)

    def lift(self, src: Path) -> int:
        """Append the alignments of ``src``; returns how many were lifted from it."""
        before = self.count
        by_name, header = self._by_name, self._header
        with pysam.AlignmentFile(str(src), check_sq=False, threads=self._threads) as bam:
            for rec in bam:
                # Edit the SAM text (RNAME, POS, RNEXT, PNEXT) and parse it back against the genome
                # header: both conversions run in htslib, ~25 % faster than to_dict/from_dict and
                # byte-identical. The record's tid cannot simply be set: pysam checks it against
                # the input (segment) header.
                f = rec.to_string().split("\t", 9)
                own = by_name.get(f[2])
                mate = own if f[6] == "=" else by_name.get(f[6])
                if own is not None:
                    f[2] = own.chrom
                    if f[3] != "0":
                        f[3] = str(int(f[3]) + own.start - 1)
                if mate is not None:
                    if f[6] != "=":
                        f[6] = mate.chrom
                    if f[7] != "0":
                        f[7] = str(int(f[7]) + mate.start - 1)
                self._out.write(pysam.AlignedSegment.fromstring("\t".join(f), header))
                self.count += 1
                if self._progress is not None and self.count % self._log_every == 0:
                    _log_lift(self._progress, self.count, self._expected)
        return self.count - before

    def close(self) -> None:
        self._out.close()

    def __enter__(self) -> AlignmentLifter:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def lift_alignments(
    sources: Sequence[Path],
    segments: Sequence[SegmentRecord],
    reference: Reference,
    dst: Path,
    *,
    threads: int = 1,
    expected: int | None = None,
    progress: str | None = None,
    log_every: int = LIFT_LOG_EVERY,
) -> int:
    """Write Mason alignments of all groups to ``dst`` (unsorted BAM) in genome coordinates."""
    with AlignmentLifter(
        segments,
        reference,
        dst,
        threads=threads,
        expected=expected,
        progress=progress,
        log_every=log_every,
    ) as lifter:
        for src in sources:
            lifter.lift(src)
    return lifter.count


def _log_lift(label: str, n: int, expected: int | None) -> None:
    if expected:
        log.info(
            "%s: lifted %s / %s alignments (%d%%)",
            label,
            f"{n:,}",
            f"{expected:,}",
            min(100, 100 * n // expected),
        )
    else:
        log.info("%s: lifted %s alignments", label, f"{n:,}")
