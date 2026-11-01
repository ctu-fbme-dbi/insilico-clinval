"""nf-core/sarek adapter (secondary analysis), built against the sarek 3.10.0 documentation.

Nextflow runs on the host (runner ``nextflow``, normally local) with ``-profile docker``; sarek
pulls its own per-process containers. All samples go into one run. Inputs:

- samplesheet ``patient,sex,status,sample,lane,fastq_1,fastq_2`` (patient = sample; sex XX/XY;
  status 0 = germline; one lane);
- ``-params-file``: ``input``, ``outdir``, ``tools``, reference (iGenomes ``genome`` or a custom
  ``fasta`` with ``igenomes_ignore``), ``skip_tools``, ``intervals``/``wes``, ``save_reference`` and
  any further ``params``;
- optional ``-c`` config with ``process.resourceLimits``.

Outputs (sarek 3.10.0 ``docs/output.md``):
``preprocessing/recalibrated/<s>/<s>.recal.cram`` (or ``markduplicates/<s>/<s>.md.cram`` when base
recalibration is skipped) and ``variant_calling/haplotypecaller/<s>/<s>.haplotypecaller.vcf.gz``
(raw) / ``….haplotypecaller.filtered.vcf.gz`` (filtered, only with filtering resources).
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from insilico_clinval.adapters import LoadedPipeline, PipelineConfigError, SampleInput, SutFailure
from insilico_clinval.schemas.enums import Sex
from insilico_clinval.schemas.pipeline import SarekConfig

PIPELINE = "nf-core/sarek"
RESERVED_PARAMS = frozenset({"input", "outdir"})
# Parameters that only say where a site keeps its reference data; ignored when adopting results.
SITE_PARAMS = frozenset({"igenomes_base", "snpeff_cache", "vep_cache"})
VERSIONS_FILE = "pipeline_info/nf_core_sarek_software_mqc_versions.yml"
DEEPVARIANT_CONFIG = "deepvariant.config"

# Sarek 3.10.0 conf/modules/deepvariant.config sets the first two arguments; the module's PAR BED
# input is always empty, so the PAR BED is passed as an argument and its directory mounted.
_DEEPVARIANT_TEMPLATE = """\
// Written by insilico_clinval: DeepVariant calls @CONTIGS@ as haploid in XY samples, except in
// the PAR regions of @PAR@.
process {
    withName: 'DEEPVARIANT_RUNDEEPVARIANT' {
        ext.args = { [
            params.wes ? "--model_type=WES" : "--model_type=WGS",
            meta.sample ? "--sample_name ${meta.sample}" : "",
            meta.sex == 'XY' ? "--haploid_contigs=@CONTIGS@ --par_regions_bed=@PAR@" : "",
        ].join(" ").trim() }
        containerOptions = "-v @PARDIR@:@PARDIR@"
    }
}
"""


@dataclass(frozen=True)
class SarekOutputs:
    alignment: Path | None
    raw_vcf: Path
    filtered_vcf: Path | None


def write_samplesheet(samples: Sequence[SampleInput], path: Path) -> Path:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"])
        for s in samples:
            sex = "XY" if s.sex == Sex.MALE else "XX"
            writer.writerow(
                [
                    s.sample_id,
                    sex,
                    0,
                    s.sample_id,
                    "lane_1",
                    s.fastq_r1.resolve(),
                    s.fastq_r2.resolve(),
                ]
            )
    return path


def sarek_params(
    cfg: SarekConfig, pipeline: LoadedPipeline, *, samplesheet: Path, outdir: Path
) -> dict[str, Any]:
    clash = sorted(RESERVED_PARAMS.intersection(cfg.params))
    if clash:
        raise PipelineConfigError(f"sarek.params must not set {clash} (set by the framework)")
    ref = cfg.reference
    params: dict[str, Any] = {
        "input": str(samplesheet.resolve()),
        "outdir": str(outdir.resolve()),
        "tools": ",".join(cfg.tools),
    }
    if cfg.skip_tools:
        params["skip_tools"] = ",".join(cfg.skip_tools)
    if ref.genome is not None:
        params["genome"] = ref.genome
    else:
        params["igenomes_ignore"] = True
    for key, param in (
        ("fasta", "fasta"), ("fasta_fai", "fasta_fai"), ("fasta_dict", "dict"), ("bwa", "bwa"),
        ("dbsnp", "dbsnp"), ("dbsnp_tbi", "dbsnp_tbi"), ("known_indels", "known_indels"),
        ("known_indels_tbi", "known_indels_tbi"),
    ):  # fmt: skip
        value = getattr(ref, key)
        if value is not None:
            params[param] = str(pipeline.resolve(value))
    if cfg.intervals is not None:
        params["intervals"] = str(pipeline.resolve(cfg.intervals))
    if cfg.wes:
        params["wes"] = True
    if cfg.save_reference:
        params["save_reference"] = True
    params.update(cfg.params)
    return params


def _groovy(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    return "'" + str(value).replace("'", "\\'") + "'"


def resource_config(limits: dict[str, Any]) -> str:
    """Nextflow config text setting ``process.resourceLimits``."""
    entries = ", ".join(f"{k}: {_groovy(v)}" for k, v in limits.items())
    return f"process {{\n    resourceLimits = [ {entries} ]\n}}\n"


def deepvariant_config(cfg: SarekConfig, pipeline: LoadedPipeline) -> str | None:
    """Nextflow config for haploid sex chromosomes in DeepVariant, or None if not configured."""
    h = cfg.deepvariant_haploid
    if h is None:
        return None
    par = pipeline.resolve(h.par_bed).resolve()
    return (
        _DEEPVARIANT_TEMPLATE.replace("@CONTIGS@", ",".join(h.contigs))
        .replace("@PARDIR@", str(par.parent))
        .replace("@PAR@", str(par))
    )


def nextflow_args(
    cfg: SarekConfig, *, params_file: Path, work_dir: Path, config_files: Sequence[Path] = ()
) -> list[str]:
    args = [
        "run", PIPELINE, "-r", cfg.release, "-profile", cfg.profile,
        "-params-file", str(params_file.resolve()), "-work-dir", str(work_dir.resolve()),
        "-ansi-log", "false",
    ]  # fmt: skip
    for config_file in config_files:
        args += ["-c", str(config_file.resolve())]
    return [*args, *cfg.extra_args]


def prepare(
    cfg: SarekConfig, pipeline: LoadedPipeline, samples: Sequence[SampleInput], run_dir: Path
) -> list[str]:
    """Write samplesheet, params file and optional config into ``run_dir``; return nextflow args."""
    run_dir.mkdir(parents=True, exist_ok=True)
    sheet = write_samplesheet(samples, run_dir / "samplesheet.csv")
    params = sarek_params(cfg, pipeline, samplesheet=sheet, outdir=run_dir / "results")
    params_file = run_dir / "params.yaml"
    params_file.write_text(yaml.safe_dump(params, sort_keys=True))
    config_files: list[Path] = []
    if cfg.resource_limits:
        config_files.append(run_dir / "resources.config")
        config_files[-1].write_text(resource_config(cfg.resource_limits))
    dv_config = deepvariant_config(cfg, pipeline)
    if dv_config is not None:
        config_files.append(run_dir / DEEPVARIANT_CONFIG)
        config_files[-1].write_text(dv_config)
    return nextflow_args(
        cfg, params_file=params_file, work_dir=run_dir / "work", config_files=config_files
    )


def collect_outputs(results: Path, sample_id: str, caller: str = "haplotypecaller") -> SarekOutputs:
    """Locate a sample's outputs in a sarek results directory."""
    calls = results / "variant_calling" / caller / sample_id
    raw = calls / f"{sample_id}.{caller}.vcf.gz"
    if not raw.is_file():
        raise SutFailure(f"sarek produced no {caller} VCF for {sample_id} (expected {raw})")
    filtered = calls / f"{sample_id}.{caller}.filtered.vcf.gz"
    alignment: Path | None = None
    for candidate in (
        results / "preprocessing" / "recalibrated" / sample_id / f"{sample_id}.recal.cram",
        results / "preprocessing" / "markduplicates" / sample_id / f"{sample_id}.md.cram",
    ):
        if candidate.is_file():
            alignment = candidate
            break
    return SarekOutputs(
        alignment=alignment, raw_vcf=raw, filtered_vcf=filtered if filtered.is_file() else None
    )


def adopt(
    cfg: SarekConfig, pipeline: LoadedPipeline, samples: Sequence[SampleInput], results: Path
) -> list[Path]:
    """Check that a sarek run made outside the framework matches ``cfg``; return its evidence.

    ``results`` is the run's ``outdir``; its parent (the launch directory) must hold the
    ``params.yaml`` and ``samplesheet.csv`` the run used (and ``deepvariant.config`` when the
    pipeline config sets haploid calling). The parameters must equal those the
    framework would write (apart from ``input``, ``outdir`` and site data locations), the release
    must match, and the samplesheet must list every sample with the simulated FASTQ files.
    """
    launch = results.parent
    params_file, sheet = launch / "params.yaml", launch / "samplesheet.csv"
    versions = results / VERSIONS_FILE
    for path in (params_file, sheet, versions):
        if not path.is_file():
            raise SutFailure(f"cannot adopt sarek results: {path} not found")
    used = yaml.safe_load(params_file.read_text()) or {}
    expected = sarek_params(cfg, pipeline, samplesheet=sheet, outdir=results)
    errors = []
    for key in sorted(set(used) | set(expected)):
        if key in RESERVED_PARAMS or key in SITE_PARAMS:
            continue
        a, b = expected.get(key), used.get(key)
        if a is None or b is None or str(a) != str(b):
            errors.append(f"{key}: pipeline config {a!r}, run {b!r}")
    release = f"{PIPELINE}: v{cfg.release.lstrip('v')}"
    if release not in versions.read_text():
        errors.append(f"release: {versions} does not record {release!r}")
    with sheet.open(newline="") as fh:
        rows = {r["sample"]: r for r in csv.DictReader(fh)}
    for s in samples:
        row = rows.get(s.sample_id)
        fastqs = (s.fastq_r1.resolve(), s.fastq_r2.resolve())
        if row is None:
            errors.append(f"sample {s.sample_id}: not in {sheet}")
        elif (Path(row["fastq_1"]).resolve(), Path(row["fastq_2"]).resolve()) != fastqs:
            errors.append(f"sample {s.sample_id}: run used other FASTQ files than the simulation")
    evidence = [params_file, sheet, versions]
    expected_dv, dv_file = deepvariant_config(cfg, pipeline), launch / DEEPVARIANT_CONFIG
    if expected_dv is not None:
        if not dv_file.is_file() or dv_file.read_text() != expected_dv:
            errors.append(f"{DEEPVARIANT_CONFIG}: missing or different from the pipeline config")
        evidence.append(dv_file)
    elif dv_file.is_file():
        errors.append(
            f"{DEEPVARIANT_CONFIG} present but the pipeline config sets no haploid calling"
        )
    if errors:
        raise SutFailure("sarek results do not match the pipeline config: " + "; ".join(errors))
    return evidence
