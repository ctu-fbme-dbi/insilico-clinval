"""Step 2 orchestration: per sample truth VCF, expected outcomes and sample spec."""

from __future__ import annotations

import logging
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

from insilico_clinval._exec import elapsed
from insilico_clinval._runner import Toolbox
from insilico_clinval.config import LoadedSuite
from insilico_clinval.genome.reference import Reference
from insilico_clinval.jsonio import write_json
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.sample import BackgroundInfo, SampleSpec
from insilico_clinval.schemas.selection import SelectionResult
from insilico_clinval.truth.background import BackgroundMerger
from insilico_clinval.truth.expected import expected_outcomes
from insilico_clinval.truth.genotype import assign_genotypes
from insilico_clinval.truth.normalize import compress_and_index, normalize_vcf
from insilico_clinval.truth.refcheck import check_ref_alleles
from insilico_clinval.truth.vcfwrite import write_truth_vcf
from insilico_clinval.variants.select import variants_per_sample

WORK_DIR = ".work"

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class SampleOutputs:
    sample_id: str
    truth_vcf: Path
    truth_tbi: Path
    expected_outcomes: Path
    sample_spec: Path
    background_removed: Path | None = None


def build_truth(
    suite: LoadedSuite,
    selection: SelectionResult,
    reference: Reference,
    out_dir: Path,
    manifest: ManifestBuilder | None,
    toolbox: Toolbox | None = None,
) -> list[SampleOutputs]:
    cfg = suite.config
    toolbox = toolbox if toolbox is not None else Toolbox()
    check_ref_alleles(selection.variants, reference)
    per_sample = variants_per_sample(selection)
    log.info("truth: %d selected variants, %d sample(s)", len(selection.variants), len(cfg.samples))
    work = out_dir / WORK_DIR
    work.mkdir(parents=True, exist_ok=True)
    fasta = suite.resolve(cfg.reference.fasta)
    outputs: list[SampleOutputs] = []
    try:
        for sample in cfg.samples:
            sid = sample.sample_id
            t_sample = time.perf_counter()
            planted = per_sample.get(sid, [])
            genotypes, warnings = assign_genotypes(cfg, sample, planted)
            merger: BackgroundMerger | None = None
            if cfg.uses_background(sample):
                assert cfg.background is not None
                log.info(
                    "%s: writing truth VCF: %d planted variants + %s background",
                    sid,
                    len(planted),
                    cfg.background.name,
                )
                merger = BackgroundMerger(
                    config=cfg.background,
                    vcf_path=suite.resolve(cfg.background.vcf),
                    bed_path=suite.resolve(cfg.background.bed),
                    reference=reference,
                    seed=cfg.suite.seed,
                    sample_id=sid,
                    planted=planted,
                    min_distance_bp=cfg.background_min_distance_bp(),
                )
            else:
                log.info("%s: writing truth VCF: %d planted variants", sid, len(planted))
            try:
                raw, n_records = write_truth_vcf(
                    work / f"{sid}.raw.vcf", sid, genotypes, reference, merger
                )
            finally:
                if merger is not None:
                    merger.close()
            log.info(
                "%s: truth VCF written in %s: %s records", sid, elapsed(t_sample), f"{n_records:,}"
            )
            background_info: BackgroundInfo | None = None
            removed_path: Path | None = None
            if merger is not None:
                assert cfg.background is not None
                removed_path = merger.write_removed(out_dir / f"{sid}.background_removed.tsv")
                outside = merger.planted_outside_bed()
                background_info = BackgroundInfo(
                    name=cfg.background.name,
                    release=cfg.background.release,
                    vcf=cfg.background.vcf,
                    bed=cfg.background.bed,
                    records=merger.written,
                    removed=len(merger.removed),
                    min_distance_bp=merger.min_distance_bp,
                    removed_log=removed_path.name,
                    planted_outside_bed=outside,
                    skipped=dict(sorted(merger.skipped.items())),
                    spanning_deletion_to_ref=merger.spanning_deletion_to_ref,
                )
                name = cfg.background.name
                warnings.append(
                    f"sample {sid}: {len(merger.removed)} {name} background records removed "
                    f"within {merger.min_distance_bp} bp of planted variants "
                    f"(see {removed_path.name})"
                )
                for reason, n in sorted(merger.skipped.items()):
                    warnings.append(f"sample {sid}: {n} {name} records skipped ({reason})")
                if outside:
                    warnings.append(
                        f"sample {sid}: {len(outside)} planted variants outside the {name} "
                        f"high-confidence BED: {', '.join(outside)}"
                    )
            log.info("%s: normalising (bcftools norm), compressing and indexing", sid)
            t0 = time.perf_counter()
            norm = normalize_vcf(
                raw,
                work / f"{sid}.norm.vcf",
                fasta,
                step=f"bcftools norm ({sid})",
                manifest=manifest,
                toolbox=toolbox,
            )
            vcf_gz = compress_and_index(norm, out_dir / f"{sid}.truth.vcf.gz")
            log.info("%s: truth VCF normalised and indexed in %s", sid, elapsed(t0))
            outcome_set, more_warnings = expected_outcomes(cfg, sample, genotypes)
            eo_path = write_json(out_dir / f"{sid}.expected_outcomes.json", outcome_set)
            spec = SampleSpec(
                schema_version=SCHEMA_VERSION,
                sample_id=sid,
                role=sample.role,
                sex=sample.sex,
                background=background_info,
                variants=sorted(
                    (g.to_sample_variant() for g in genotypes), key=lambda s: s.variant_id
                ),
                truth_vcf=vcf_gz.name,
                expected_outcomes=eo_path.name,
            )
            spec_path = write_json(out_dir / f"{sid}.sample.json", spec)
            log.info("%s: truth finished in %s", sid, elapsed(t_sample))
            if manifest is not None:
                for w in warnings + more_warnings:
                    manifest.warn(w)
            outputs.append(
                SampleOutputs(
                    sample_id=sid,
                    truth_vcf=vcf_gz,
                    truth_tbi=Path(f"{vcf_gz}.tbi"),
                    expected_outcomes=eo_path,
                    sample_spec=spec_path,
                    background_removed=removed_path,
                )
            )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return outputs
