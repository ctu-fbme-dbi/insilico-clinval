"""Step 3b orchestration for one sample: Mason per group or batch, merge, sort, check."""

from __future__ import annotations

import logging
import time
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from insilico_clinval._exec import elapsed
from insilico_clinval._runner import Toolbox
from insilico_clinval.genome.reference import Reference
from insilico_clinval.jsonio import write_json
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import PlantingStatus, SimulationMode
from insilico_clinval.schemas.sample import SampleSpec
from insilico_clinval.schemas.simulation import MasonRun, SimulationProfile, SimulationSpec
from insilico_clinval.schemas.suite import SimulationConfig
from insilico_clinval.simprep.prepare import SimPrep, read_planted
from insilico_clinval.simulate import SimulationError
from insilico_clinval.simulate.check import planting_check, write_check_tsv
from insilico_clinval.simulate.mason import fragment_count, mason_args, mason_seed
from insilico_clinval.simulate.merge import AlignmentLifter, FastqGzConcat

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class SimulationOutputs:
    sample_id: str
    fastq_r1: Path
    fastq_r2: Path
    true_bam: Path
    true_bai: Path
    spec: Path
    check_json: Path
    check_tsv: Path
    complete: bool


def sim_dir_name(sample_id: str) -> str:
    return f"{sample_id}.sim"


def simulate_sample(
    *,
    prep: SimPrep,
    spec: SampleSpec,
    sim: SimulationConfig,
    profile: SimulationProfile,
    profile_sha256: str,
    seed: int,
    reference: Reference,
    out_dir: Path,
    toolbox: Toolbox,
    manifest: ManifestBuilder | None,
    keep_work: bool = False,
) -> SimulationOutputs:
    sid = prep.sample_id
    sim_dir = prep.sim_dir.resolve()
    out_dir = out_dir.resolve()
    runs: list[MasonRun] = []
    t_sample = time.perf_counter()
    fragments_of = {g.name: fragment_count(profile, g.length, g.ploidy) for g in prep.groups}
    expected = 2 * sum(fragments_of.values())
    log.info(
        "%s: simulating %d Mason run(s), %gx, 2 x %d bp, %d Mason thread(s), %s alignments",
        sid,
        len(prep.groups),
        profile.coverage,
        profile.read_length,
        profile.threads,
        f"{expected:,}",
    )
    # Compression, lifting and sorting use the profile's threads too; none of these outputs
    # depends on the thread count.
    threads = profile.threads
    fastq_r1 = out_dir / f"{sid}.R1.fastq.gz"
    fastq_r2 = out_dir / f"{sid}.R2.fastq.gz"
    unsorted = sim_dir / "true.unsorted.bam"
    # Each Mason run is merged into the outputs (FASTQ appended, alignments lifted) and its raw
    # files deleted as soon as it ends, so the disk peak stays near that of the largest run (see
    # simulation.mason_batch_bp). The merge of run k runs in one background thread, in run order,
    # while Mason simulates run k+1; the outputs are the same as merging one run after another.
    with (
        FastqGzConcat(fastq_r1, threads=threads) as gz1,
        FastqGzConcat(fastq_r2, threads=threads) as gz2,
        AlignmentLifter(
            prep.segments, reference, unsorted, threads=threads, expected=expected, progress=sid
        ) as lifter,
        ThreadPoolExecutor(max_workers=1, thread_name_prefix="merge") as merger,
    ):

        def merge_run(name: str, r1: Path, r2: Path, bam: Path) -> None:
            t0 = time.perf_counter()
            gz1.append(r1)
            gz2.append(r2)
            log.info("%s: %s: FASTQ appended in %s", sid, name, elapsed(t0))
            t0 = time.perf_counter()
            n = lifter.lift(bam)
            log.info(
                "%s: %s: %s true alignments lifted in %s (%s / %s)",
                sid,
                name,
                f"{n:,}",
                elapsed(t0),
                f"{lifter.count:,}",
                f"{expected:,}",
            )
            if not keep_work:
                for raw in (r1, r2, bam):
                    raw.unlink(missing_ok=True)

        merges: list[Future[None]] = []
        for group in prep.groups:
            for done in merges:  # fail early if an earlier merge failed
                if done.done():
                    done.result()
            name = group.name
            fragments = fragments_of[name]
            group_seed = mason_seed(seed, sid, name)
            r1, r2, bam = f"{name}_1.fq", f"{name}_2.fq", f"{name}.true.bam"
            args = mason_args(
                profile,
                fasta=group.fasta.name,
                vcf=group.vcf.name if group.vcf is not None else None,
                fragments=fragments,
                seed=group_seed,
                read_name_prefix=f"{sid}.{name}.",
                out_r1=r1,
                out_r2=r2,
                out_alignment=bam,
            )
            log.info(
                "%s: mason_simulator %s started: %s fragments over %s bp",
                sid,
                name,
                f"{fragments:,}",
                f"{group.length:,}",
            )
            t0 = time.perf_counter()
            # Mason's stderr goes straight to its log: followable with tail -f, kept if it fails.
            mason_log = sim_dir / f"mason_{name}.log"
            log.info("%s: Mason output: %s", sid, mason_log)
            toolbox.run(
                "mason_simulator",
                args,
                step=f"mason_simulator ({sid}, {name})",
                manifest=manifest,
                mounts=[sim_dir],
                cwd=sim_dir,
                stderr_path=mason_log,
                # Mason keeps all reads in (unlinked) temporary files until it joins them at the
                # end; by default they go to /tmp, i.e. the container's (or system) disk. A full
                # genome filled it and Mason died with SIGSEGV (it does not check write errors).
                env={"TMPDIR": str(sim_dir)},
            )
            log.info("%s: mason_simulator %s done in %s", sid, name, elapsed(t0))
            if sim.mode == SimulationMode.GENOME and not keep_work:
                group.fasta.unlink(missing_ok=True)
                Path(f"{group.fasta}.fai").unlink(missing_ok=True)

            merges.append(merger.submit(merge_run, name, sim_dir / r1, sim_dir / r2, sim_dir / bam))
            runs.append(
                MasonRun(
                    group=group.group,
                    name=name if group.batched else None,
                    segment_names=[s.name for s in group.segments] if group.batched else None,
                    ploidy=group.ploidy,
                    segments=len(group.segments),
                    length=group.length,
                    fragments=fragments,
                    seed=group_seed,
                    vcf=_rel(group.vcf, out_dir) if group.vcf is not None else None,
                    vcf_records=group.records,
                    args=args,
                )
            )
        if merges and not merges[-1].done():
            log.info("%s: waiting for the last merge (FASTQ + lift)", sid)
        for future in merges:
            future.result()

    true_bam = out_dir / f"{sid}.true.bam"
    log.info("%s: sorting and indexing the true BAM (samtools)", sid)
    t0 = time.perf_counter()
    toolbox.run(
        "samtools",
        ["sort", "--no-PG", "-@", str(threads), "-o", str(true_bam), str(unsorted)],
        step=f"samtools sort ({sid})",
        manifest=manifest,
        mounts=[out_dir],
    )
    if not keep_work:
        unsorted.unlink(missing_ok=True)
    toolbox.run(
        "samtools",
        ["index", "-@", str(threads), str(true_bam)],
        step=f"samtools index ({sid})",
        manifest=manifest,
        mounts=[out_dir],
    )
    true_bai = Path(f"{true_bam}.bai")
    log.info("%s: true BAM sorted and indexed in %s", sid, elapsed(t0))

    check = planting_check(
        true_bam, reference, sid, prep.planted, spec.variants, min_alt_reads=sim.min_alt_reads
    )
    n_ok = sum(v.status == PlantingStatus.OK for v in check.variants)
    log.info("%s: planting check %d/%d variants ok", sid, n_ok, len(check.variants))
    check_json = write_json(out_dir / f"{sid}.planting_check.json", check)
    check_tsv = write_check_tsv(check, out_dir / f"{sid}.planting_check.tsv")

    sim_spec = SimulationSpec(
        schema_version=SCHEMA_VERSION,
        sample_id=sid,
        mode=sim.mode,
        profile=profile,
        profile_sha256=profile_sha256,
        flank_bp=sim.flank_bp if sim.mode == SimulationMode.REGION else None,
        segments=_rel(prep.segments_tsv, out_dir),
        runs=runs,
        decomposed=prep.decomposed,
        background_dropped=len(prep.dropped),
        background_dropped_log=_rel(prep.dropped_log, out_dir) if prep.dropped_log else None,
        fastq_r1=fastq_r1.name,
        fastq_r2=fastq_r2.name,
        true_bam=true_bam.name,
        planting_check=check_json.name,
        planting_complete=check.complete,
    )
    spec_path = write_json(out_dir / f"{sid}.simulation.json", sim_spec)
    log.info("%s: simulation finished in %s", sid, elapsed(t_sample))
    return SimulationOutputs(
        sample_id=sid,
        fastq_r1=fastq_r1,
        fastq_r2=fastq_r2,
        true_bam=true_bam,
        true_bai=true_bai,
        spec=spec_path,
        check_json=check_json,
        check_tsv=check_tsv,
        complete=check.complete,
    )


@dataclass(frozen=True)
class RecheckOutputs:
    check_json: Path
    check_tsv: Path
    spec: Path
    complete: bool


def recheck_planting(
    out_dir: Path,
    spec: SampleSpec,
    truth_vcf: Path,
    reference: Reference,
    *,
    min_alt_reads: int,
) -> RecheckOutputs:
    """Repeat the planting check of a finished simulation from its true BAM (no Mason run).

    Rewrites ``<sample>.planting_check.json``/``.tsv`` and ``planting_complete`` in
    ``<sample>.simulation.json``; the reads and the true BAM are not touched.
    """
    sid = spec.sample_id
    spec_path = out_dir / f"{sid}.simulation.json"
    if not spec_path.is_file():
        raise SimulationError(f"{spec_path} not found (run `clinval simulate` first)")
    sim_spec = SimulationSpec.model_validate_json(spec_path.read_text())
    true_bam = out_dir / sim_spec.true_bam
    if not true_bam.is_file() or not Path(f"{true_bam}.bai").is_file():
        raise SimulationError(f"{true_bam} or its index not found")
    planted = read_planted(truth_vcf, [v.variant_id for v in spec.variants], reference)
    check = planting_check(
        true_bam, reference, sid, planted, spec.variants, min_alt_reads=min_alt_reads
    )
    n_ok = sum(v.status == PlantingStatus.OK for v in check.variants)
    log.info("%s: planting check (recheck) %d/%d variants ok", sid, n_ok, len(check.variants))
    check_json = write_json(out_dir / f"{sid}.planting_check.json", check)
    check_tsv = write_check_tsv(check, out_dir / f"{sid}.planting_check.tsv")
    updated = sim_spec.model_copy(
        update={"planting_check": check_json.name, "planting_complete": check.complete}
    )
    write_json(spec_path, updated)
    return RecheckOutputs(check_json, check_tsv, spec_path, check.complete)


def _rel(path: Path, base: Path) -> str:
    return str(path.resolve().relative_to(base.resolve()))
