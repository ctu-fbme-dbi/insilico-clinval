"""Step 3a orchestration: segments, segment FASTA, Mason VCF per ploidy group.

The truth VCF is read twice: first the planted variants (located through their ids, which carry
chrom and position, and tabix; they define region-mode segments and are reserved on their
haplotypes first), then the records to write — fetched per segment in region mode, streamed in
genome mode — so background genomes of millions of records are never held in memory. Planted
variants that cannot be represented are errors; background records that cannot be represented
are dropped and logged.
"""

from __future__ import annotations

import heapq
import logging
import time
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pysam

from insilico_clinval._exec import elapsed
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.enums import Sex, SimulationMode
from insilico_clinval.schemas.simulation import SegmentRecord
from insilico_clinval.schemas.suite import Region, SimulationConfig
from insilico_clinval.simprep import SimPrepError
from insilico_clinval.simprep.events import HaplotypeEvents, Part, decompose
from insilico_clinval.simprep.mason_vcf import MasonVcfWriter, validate_mason_vcf
from insilico_clinval.simprep.segments import (
    SegmentIndex,
    build_segments,
    genome_intervals,
    read_bed_intervals,
    region_intervals,
)

log = logging.getLogger(__name__)

Group = Literal["diploid", "haploid"]
FASTA_LINE = 60


def group_of(ploidy: int) -> Group:
    return "haploid" if ploidy == 1 else "diploid"


@dataclass(frozen=True)
class TruthRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple[int, ...]
    variant_id: str | None  # None for background records

    @property
    def end(self) -> int:
        return self.pos + len(self.ref) - 1


@dataclass(frozen=True)
class GroupInputs:
    """Inputs of one Mason run: a ploidy group, or one batch of it."""

    group: Group
    ploidy: Literal[1, 2]
    segments: list[SegmentRecord]
    fasta: Path
    vcf: Path | None
    records: int
    name: str  # file prefix, seed and read-name label: the group, or `<group>-<n>` for a batch
    batched: bool = False

    @property
    def length(self) -> int:
        return sum(s.end - s.start + 1 for s in self.segments)


@dataclass
class SimPrep:
    sample_id: str
    sim_dir: Path
    segments: list[SegmentRecord]
    segments_tsv: Path
    groups: list[GroupInputs]
    planted: list[TruthRecord]
    decomposed: list[str]
    dropped: list[tuple[TruthRecord, str]] = field(default_factory=list)
    dropped_log: Path | None = None


def batch_segments(segments: list[SegmentRecord], max_bp: int | None) -> list[list[SegmentRecord]]:
    """Consecutive segments packed greedily into batches of at most ``max_bp`` bp.

    A segment longer than ``max_bp`` forms a batch of its own; ``None`` gives one batch.
    """
    if max_bp is None:
        return [segments]
    batches: list[list[SegmentRecord]] = []
    current: list[SegmentRecord] = []
    current_bp = 0
    for seg in segments:
        n = seg.end - seg.start + 1
        if current and current_bp + n > max_bp:
            batches.append(current)
            current, current_bp = [], 0
        current.append(seg)
        current_bp += n
    if current:
        batches.append(current)
    return batches


def _truth_record(r: pysam.VariantRecord, sample: str, path: Path) -> TruthRecord | None:
    """A biallelic truth record with its GT; None if the genotype is missing."""
    alts = r.alts or ()
    if len(alts) != 1 or r.ref is None:
        raise SimPrepError(f"{path.name}: record {r.chrom}:{r.pos} is not biallelic")
    gt = r.samples[sample]["GT"]
    if gt is None or any(a is None for a in gt):
        return None
    return TruthRecord(
        chrom=r.chrom,
        pos=r.pos,
        ref=r.ref.upper(),
        alt=alts[0].upper(),
        gt=tuple(int(a) for a in gt),
        variant_id=None if "TV_BG" in r.info else str(r.info["TV_ID"]),
    )


def read_truth(path: Path, *, planted_only: bool = False) -> list[TruthRecord]:
    return [r for r in iter_truth(path) if not planted_only or r.variant_id is not None]


def iter_truth(path: Path, segments: list[SegmentRecord] | None = None) -> Iterator[TruthRecord]:
    """Truth records in file order, or (with ``segments``) those overlapping the segments.

    Segments must be sorted in FASTA order; a record overlapping two segments is yielded once.
    """
    with pysam.VariantFile(str(path)) as vcf:
        sample = next(iter(vcf.header.samples))
        if segments is None:
            records: Iterator[pysam.VariantRecord] = iter(vcf)
        else:
            records = _fetch_segments(vcf, segments)
        for r in records:
            record = _truth_record(r, sample, path)
            if record is not None:
                yield record


def _fetch_segments(
    vcf: pysam.VariantFile, segments: list[SegmentRecord]
) -> Iterator[pysam.VariantRecord]:
    contigs = set(vcf.header.contigs)
    seen_crossing: set[tuple[str, int, str, str]] = set()
    for seg in segments:
        if seg.chrom not in contigs:
            continue
        for r in vcf.fetch(seg.chrom, seg.start - 1, seg.end):
            if r.pos < seg.start:
                key = (r.chrom, r.pos, str(r.ref), str(r.alts))
                if key in seen_crossing:
                    continue
                seen_crossing.add(key)
            yield r


def read_planted(path: Path, variant_ids: list[str], reference: Reference) -> list[TruthRecord]:
    """Planted records located by their ids (``{category}:{chrom}:{pos}:{ref}:{alt}``) + tabix."""
    found: list[TruthRecord] = []
    missing: list[str] = []
    with pysam.VariantFile(str(path)) as vcf:
        sample = next(iter(vcf.header.samples))
        for vid in variant_ids:
            _, chrom, pos_text, _, _ = vid.split(":")
            pos = int(pos_text)
            hits = [
                _truth_record(r, sample, path)
                for r in vcf.fetch(chrom, pos - 1, pos)
                if r.pos == pos and "TV_ID" in r.info and r.info["TV_ID"] == vid
            ]
            record = hits[0] if hits else None
            if record is None:
                missing.append(vid)
            else:
                found.append(record)
    if missing:
        raise SimPrepError(f"{path.name}: planted variants not found: {', '.join(missing)}")
    return sorted(found, key=lambda r: reference.sort_key(r.chrom, r.pos))


def _intervals(
    sim: SimulationConfig,
    planted: list[TruthRecord],
    reference: Reference,
    extra_bed: Path | None,
) -> list[tuple[str, int, int]]:
    if sim.mode == SimulationMode.GENOME:
        return genome_intervals(reference, sim.genome_contigs)
    extra = read_bed_intervals(extra_bed, reference) if extra_bed is not None else []
    return region_intervals(
        ((r.chrom, r.pos, r.ref) for r in planted), reference, sim.flank_bp, extra
    )


def _write_fasta(path: Path, segments: list[SegmentRecord], reference: Reference) -> None:
    with path.open("w") as fh:
        for seg in segments:
            seq = reference.fetch(seg.chrom, seg.start - 1, seg.end)
            fh.write(f">{seg.name}\n")
            for i in range(0, len(seq), FASTA_LINE):
                fh.write(seq[i : i + FASTA_LINE] + "\n")
    pysam.faidx(str(path))


def _ploidy_gt(record: TruthRecord, ploidy: int) -> tuple[int, ...] | str:
    """GT adapted to the segment ploidy, or the reason why it cannot be represented."""
    if ploidy == 2:
        return record.gt if len(record.gt) == 2 else "haploid_gt_on_diploid_segment"
    if len(record.gt) == 1:
        return record.gt
    if len(set(record.gt)) == 1:
        return record.gt[:1]
    return "het_on_haploid_segment"


def _indel_at_segment_end(record: TruthRecord, seg: SegmentRecord) -> bool:
    """An indel whose REF span reaches the last base of the segment.

    mason_simulator (seqan-v2.4.0) corrupts its heap on such a record (std::bad_alloc, exit 139),
    seen in region mode, where segments end inside the chromosome.
    """
    return len(record.ref) != len(record.alt) and record.end >= seg.end


def _haplotypes(gt: tuple[int, ...]) -> list[int]:
    """0-based haplotype indices carrying the ALT allele."""
    return [i for i, allele in enumerate(gt) if allele != 0]


class _SortedParts:
    """Buffers parts so each group's Mason VCF stays sorted (decomposed parts may interleave).

    At one POS an SNV is written before an indel: a delins with a different first base
    decomposes into an SNV at POS plus an indel anchored at POS, and the indel changes only the
    bases after its anchor, so this is the order in which the parts' spans start.
    """

    def __init__(self, writer: MasonVcfWriter, order: dict[str, int]):
        self.writer = writer
        self.order = order
        self.heap: list[tuple[int, int, bool, str, str, str, tuple[int, ...], str]] = []

    def push(self, seg: SegmentRecord, part: Part, gt: tuple[int, ...], vid: str) -> None:
        pos = part.pos - seg.start + 1
        item = (self.order[seg.name], pos, part.is_indel, part.ref, part.alt, seg.name, gt, vid)
        heapq.heappush(self.heap, item)

    def flush_before(self, seg_order: int, pos: int) -> None:
        while self.heap and (self.heap[0][0], self.heap[0][1]) < (seg_order, pos):
            self._pop()

    def flush(self) -> None:
        while self.heap:
            self._pop()

    def _pop(self) -> None:
        _, pos, _, ref, alt, name, gt, vid = heapq.heappop(self.heap)
        self.writer.write(name, pos, ref, alt, gt, vid)


def prepare_sample(
    *,
    sample_id: str,
    sex: Sex,
    truth_vcf: Path,
    planted_ids: list[str],
    reference: Reference,
    pars: tuple[Region, ...],
    sim: SimulationConfig,
    extra_bed: Path | None,
    min_segment_bp: int,
    sim_dir: Path,
) -> SimPrep:
    """Write Mason inputs of one sample into ``sim_dir``."""
    sim_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    log.info("%s: preparing Mason inputs (%s mode)", sample_id, sim.mode.value)
    planted = read_planted(truth_vcf, planted_ids, reference)
    segments = build_segments(_intervals(sim, planted, reference, extra_bed), reference, sex, pars)
    short = [s.name for s in segments if s.end - s.start + 1 < min_segment_bp]
    if short:
        raise SimPrepError(
            f"segments shorter than the longest fragment ({min_segment_bp} bp): "
            f"{', '.join(short)} (increase simulation.flank_bp)"
        )
    index = SegmentIndex(segments)
    events: dict[tuple[str, int], HaplotypeEvents] = {}
    decomposed: list[str] = []
    errors: list[str] = []
    for r in planted:
        seg = index.overlapping(r.chrom, r.pos, r.end)
        if seg is None or not index.contains(r.chrom, r.pos, r.end):
            errors.append(f"{r.variant_id}: not inside a simulated segment")
            continue
        gt = _ploidy_gt(r, seg.ploidy)
        if isinstance(gt, str):
            errors.append(f"{r.variant_id}: GT {r.gt} on a {seg.ploidy}-ploid segment ({gt})")
            continue
        if _indel_at_segment_end(r, seg):
            errors.append(f"{r.variant_id}: indel reaches the last base of its segment")
            continue
        parts = decompose(r.pos, r.ref, r.alt)
        if parts != [Part(r.pos, r.ref, r.alt)]:
            decomposed.append(str(r.variant_id))
        for hap in _haplotypes(gt):
            hap_events = events.setdefault((seg.name, hap), HaplotypeEvents(planted=[]))
            for part in parts:
                other = hap_events.conflict_with(part)
                if other is not None:
                    errors.append(
                        f"{r.variant_id}: overlaps another planted variant on haplotype "
                        f"{hap + 1} ({other.pos} {other.ref}>{other.alt})"
                    )
                hap_events.planted.append(part)
    if errors:
        raise SimPrepError("planted variants cannot be simulated:\n  " + "\n  ".join(errors))

    groups: list[GroupInputs] = []
    order = {s.name: i for i, s in enumerate(segments)}
    group_segments: dict[int, list[SegmentRecord]] = {}
    for seg in segments:
        group_segments.setdefault(seg.ploidy, []).append(seg)
    # One Mason run per ploidy group, or per batch of it (simulation.mason_batch_bp).
    runs: list[tuple[str, int, list[SegmentRecord]]] = []
    for ploidy in sorted(group_segments, reverse=True):
        batches = batch_segments(group_segments[ploidy], sim.mason_batch_bp)
        for i, batch in enumerate(batches, start=1):
            name = group_of(ploidy) if sim.mason_batch_bp is None else f"{group_of(ploidy)}-{i}"
            runs.append((name, ploidy, batch))
    run_of = {seg.name: name for name, _, batch in runs for seg in batch}
    writers: dict[str, _SortedParts] = {}
    counts: dict[str, int] = {}
    for name, _, batch in runs:
        _write_fasta(sim_dir / f"{name}.fa", batch, reference)
        writer = MasonVcfWriter(sim_dir / f"{name}.mason.vcf", sample_id, batch)
        writers[name] = _SortedParts(writer, order)
        counts[name] = 0

    dropped: list[tuple[TruthRecord, str]] = []
    try:
        region = sim.mode == SimulationMode.REGION
        for r in iter_truth(truth_vcf, segments if region else None):
            seg = index.overlapping(r.chrom, r.pos, r.end)
            if seg is None:
                continue  # region mode: outside the simulated regions
            run = run_of[seg.name]
            sink = writers[run]
            sink.flush_before(order[seg.name], r.pos - seg.start + 1)
            if r.variant_id is not None:
                gt = _ploidy_gt(r, seg.ploidy)
                assert not isinstance(gt, str)
                for part in decompose(r.pos, r.ref, r.alt):
                    sink.push(seg, part, gt, r.variant_id)
                    counts[run] += 1
                continue
            if not index.contains(r.chrom, r.pos, r.end):
                dropped.append((r, "crosses_segment_boundary"))
                continue
            if _indel_at_segment_end(r, seg):
                dropped.append((r, "indel_at_segment_end"))
                continue
            gt_or_reason = _ploidy_gt(r, seg.ploidy)
            if isinstance(gt_or_reason, str):
                dropped.append((r, gt_or_reason))
                continue
            parts = decompose(r.pos, r.ref, r.alt)
            haps = _haplotypes(gt_or_reason)
            carriers = [events.setdefault((seg.name, h), HaplotypeEvents(planted=[])) for h in haps]
            if any(e.background_conflict(p) for e in carriers for p in parts):
                dropped.append((r, "overlap_on_haplotype"))
                continue
            for e in carriers:
                for p in parts:
                    e.accept_background(p)
            for part in parts:
                sink.push(seg, part, gt_or_reason, ".")
                counts[run] += 1
        for sink in writers.values():
            sink.flush()
    finally:
        for sink in writers.values():
            sink.writer.close()

    for name, ploidy, batch in runs:
        vcf_path: Path | None = sim_dir / f"{name}.mason.vcf"
        assert vcf_path is not None
        if counts[name] == 0:
            vcf_path.unlink()
            vcf_path = None
        else:
            validate_mason_vcf(vcf_path, sim_dir / f"{name}.fa")
        groups.append(
            GroupInputs(
                group=group_of(ploidy),
                ploidy=1 if ploidy == 1 else 2,
                segments=batch,
                fasta=sim_dir / f"{name}.fa",
                vcf=vcf_path,
                records=counts[name],
                name=name,
                batched=sim.mason_batch_bp is not None,
            )
        )

    segments_tsv = sim_dir / "segments.tsv"
    with segments_tsv.open("w") as fh:
        fh.write("segment\tchrom\tstart\tend\tploidy\tgroup\n")
        for seg in segments:
            fh.write(
                f"{seg.name}\t{seg.chrom}\t{seg.start}\t{seg.end}\t{seg.ploidy}\t"
                f"{group_of(seg.ploidy)}\n"
            )
    dropped_log: Path | None = None
    if dropped:
        dropped_log = sim_dir / "background_dropped.tsv"
        with dropped_log.open("w") as fh:
            fh.write("chrom\tpos\tref\talt\tgt\treason\n")
            for r, reason in dropped:
                gt_text = "|".join(str(a) for a in r.gt)
                fh.write(f"{r.chrom}\t{r.pos}\t{r.ref}\t{r.alt}\t{gt_text}\t{reason}\n")
    log.info(
        "%s: Mason inputs ready in %s: %d segments, %s records (%s), %d background dropped",
        sample_id,
        elapsed(t0),
        len(segments),
        f"{sum(g.records for g in groups):,}",
        ", ".join(f"{g.name} {g.records:,}" for g in groups),
        len(dropped),
    )
    return SimPrep(
        sample_id=sample_id,
        sim_dir=sim_dir,
        segments=segments,
        segments_tsv=segments_tsv,
        groups=groups,
        planted=planted,
        decomposed=sorted(decomposed),
        dropped=dropped,
        dropped_log=dropped_log,
    )
