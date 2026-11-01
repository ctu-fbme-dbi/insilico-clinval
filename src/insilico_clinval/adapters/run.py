"""Run a pipeline under test on the simulated samples of a run directory.

Layout (``sut_dir`` = ``<run>/sut/<pipeline name>/``)::

    secondary/                 sarek: samplesheet.csv, params.yaml, results/, work/
    secondary/<sample>/        command adapter work directory
    tertiary/<sample>/         tertiary input VCF, annotator/classifier outputs
    <sample>.secondary.json    SecondaryResult
    <sample>.findings.json     FindingSet
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from insilico_clinval._runner import Toolbox
from insilico_clinval.adapters import LoadedPipeline, SampleInput, SutFailure
from insilico_clinval.adapters import sarek as sarek_adapter
from insilico_clinval.adapters.command import run_black_box, substitute
from insilico_clinval.adapters.tertiary import (
    prepare_input,
    run_bias2015,
    run_command_tertiary,
)
from insilico_clinval.jsonio import write_json
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.finding import FindingSet, SecondaryResult


def sut_dir(run_dir: Path, pipeline: LoadedPipeline) -> Path:
    return run_dir / "sut" / pipeline.config.name


def _rel(path: Path | None, base: Path) -> str | None:
    """Path relative to the run directory, or absolute when outside it."""
    if path is None:
        return None
    resolved, root = path.resolve(), base.resolve()
    return str(resolved.relative_to(root)) if resolved.is_relative_to(root) else str(resolved)


def _abs(value: str | None, base: Path) -> Path | None:
    if value is None:
        return None
    p = Path(value)
    return p if p.is_absolute() else base / p


def run_secondary(
    pipeline: LoadedPipeline,
    samples: Sequence[SampleInput],
    *,
    run_dir: Path,
    fasta: Path,
    toolbox: Toolbox,
    manifest: ManifestBuilder | None,
    sarek_results: Path | None = None,
) -> list[tuple[SecondaryResult, Path]]:
    """Run the secondary stage; with ``sarek_results``, adopt a sarek run made outside."""
    cfg = pipeline.config.secondary
    out = sut_dir(run_dir, pipeline)
    work = out / "secondary"
    work.mkdir(parents=True, exist_ok=True)
    name = pipeline.config.name
    found: dict[str, tuple[str, str, Path | None, Path, Path | None]] = {}

    if sarek_results is not None and cfg.adapter != "sarek":
        raise SutFailure(f"sarek results given, but the secondary adapter is {cfg.adapter}")
    if cfg.adapter == "sarek":
        assert cfg.sarek is not None
        if manifest is not None:
            manifest.add_tool(sarek_adapter.PIPELINE, cfg.sarek.release, None)
        if sarek_results is None:
            args = sarek_adapter.prepare(cfg.sarek, pipeline, samples, work)
            step = f"nf-core/sarek ({name})"
            toolbox.run("nextflow", args, step=step, manifest=manifest, cwd=work)
            results_dir = work / "results"
        else:
            evidence = sarek_adapter.adopt(cfg.sarek, pipeline, samples, sarek_results)
            results_dir = sarek_results
            if manifest is not None:
                for path in evidence:
                    manifest.add_input(f"sarek_run:{path.name}", path)
                manifest.warn(f"{name}: sarek results adopted from {sarek_results}")
        for s in samples:
            o = sarek_adapter.collect_outputs(results_dir, s.sample_id, cfg.sarek.variant_caller)
            found[s.sample_id] = (
                sarek_adapter.PIPELINE, cfg.sarek.release, o.alignment, o.raw_vcf, o.filtered_vcf
            )  # fmt: skip
    elif cfg.adapter == "command":
        assert cfg.command is not None
        spec = cfg.command
        for s in samples:
            sample_dir = work / s.sample_id
            values = {
                "sample": s.sample_id,
                "fastq_r1": str(s.fastq_r1.resolve()),
                "fastq_r2": str(s.fastq_r2.resolve()),
                "fasta": str(fasta.resolve()),
                "outdir": str(sample_dir.resolve()),
            }
            run_black_box(
                toolbox,
                spec.command,
                values,
                workdir=sample_dir,
                step=f"{spec.command.tool_name} ({s.sample_id})",
                manifest=manifest,
                mounts=[s.fastq_r1.parent, fasta.parent],
            )

            def output(
                template: str | None, values: dict[str, str] = values, base: Path = sample_dir
            ) -> Path | None:
                return None if template is None else base / substitute(template, values)

            raw = output(spec.outputs.raw_vcf)
            assert raw is not None
            if not raw.is_file():
                raise SutFailure(f"{spec.command.tool_name} wrote no raw VCF {raw}")
            alignment = output(spec.outputs.alignment)
            filtered = output(spec.outputs.filtered_vcf)
            found[s.sample_id] = (
                spec.command.tool_name,
                spec.command.version,
                alignment if alignment is not None and alignment.is_file() else None,
                raw,
                filtered if filtered is not None and filtered.is_file() else None,
            )
    else:
        assert cfg.precomputed is not None
        pre = cfg.precomputed
        for s in samples:
            files = pre.samples.get(s.sample_id)
            if files is None:
                raise SutFailure(f"precomputed outputs missing for sample {s.sample_id}")
            raw = pipeline.resolve(files.raw_vcf)
            if not raw.is_file():
                raise SutFailure(f"precomputed raw VCF not found: {raw}")
            found[s.sample_id] = (
                pre.tool_name,
                pre.version,
                pipeline.resolve(files.alignment) if files.alignment else None,
                raw,
                pipeline.resolve(files.filtered_vcf) if files.filtered_vcf else None,
            )
        if manifest is not None:
            manifest.add_tool(pre.tool_name, pre.version, None)

    results: list[tuple[SecondaryResult, Path]] = []
    for s in samples:
        tool, version, alignment, raw, filtered = found[s.sample_id]
        result = SecondaryResult(
            schema_version=SCHEMA_VERSION,
            sample_id=s.sample_id,
            pipeline=name,
            adapter=cfg.adapter,
            tool=tool,
            tool_version=version,
            alignment=_rel(alignment, run_dir),
            raw_vcf=_rel(raw, run_dir) or "",
            filtered_vcf=_rel(filtered, run_dir),
        )
        path = write_json(out / f"{s.sample_id}.secondary.json", result)
        results.append((result, path))
    return results


def load_secondary(run_dir: Path, pipeline: LoadedPipeline, sample_id: str) -> SecondaryResult:
    path = sut_dir(run_dir, pipeline) / f"{sample_id}.secondary.json"
    if not path.is_file():
        raise SutFailure(f"{path} not found (run the secondary stage first)")
    return SecondaryResult.model_validate_json(path.read_text())


def run_tertiary(
    pipeline: LoadedPipeline,
    secondary: Sequence[SecondaryResult],
    *,
    run_dir: Path,
    fasta: Path,
    toolbox: Toolbox,
    manifest: ManifestBuilder | None,
) -> list[tuple[FindingSet, Path]]:
    cfg = pipeline.config.tertiary
    out = sut_dir(run_dir, pipeline)
    results: list[tuple[FindingSet, Path]] = []
    if cfg.adapter == "none":
        return results
    for sec in secondary:
        sid = sec.sample_id
        if cfg.input == "filtered" and sec.filtered_vcf is not None:
            source = _abs(sec.filtered_vcf, run_dir)
        else:
            source = _abs(sec.raw_vcf, run_dir)
            if cfg.input == "filtered" and manifest is not None:
                manifest.warn(f"{pipeline.config.name}/{sid}: no filtered VCF; raw VCF used")
        assert source is not None
        workdir = out / "tertiary" / sid
        vcf = prepare_input(
            toolbox, cfg, pipeline, source, workdir, sample_id=sid, manifest=manifest
        )
        if cfg.adapter == "bias2015":
            findings = run_bias2015(
                toolbox, cfg, pipeline, vcf, workdir, sample_id=sid, manifest=manifest
            )
            tool, version = "BIAS-2015", cfg.bias2015.version if cfg.bias2015 else ""
        else:
            assert cfg.command is not None
            findings = run_command_tertiary(
                toolbox, cfg, vcf, workdir, sample_id=sid, fasta=fasta, manifest=manifest
            )
            tool, version = cfg.command.command.tool_name, cfg.command.command.version
        finding_set = FindingSet(
            schema_version=SCHEMA_VERSION,
            sample_id=sid,
            pipeline=pipeline.config.name,
            adapter=cfg.adapter,
            tool=tool,
            tool_version=version,
            input_vcf=_rel(vcf, run_dir) or "",
            findings=sorted(findings, key=lambda f: (f.chrom, f.pos, f.ref, f.alt)),
        )
        path = write_json(out / f"{sid}.findings.json", finding_set)
        results.append((finding_set, path))
    return results
