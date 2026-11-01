"""Tertiary analysis: input preparation, BIAS-2015 (via Nirvana or VEP) and black-box classifiers.

Input preparation is part of the system under test (the laboratory's policy): the filtered or raw
VCF, optionally PASS records only and only the reportable regions (``bcftools view -f PASS -T``),
split into biallelic records (``bcftools norm -m -any``).

BIAS-2015 (v3.0.0) runs in three tools: an annotator writes the JSON BIAS reads (``nirvana``, or
``vep`` with the custom annotations and plugins of BIAS's doc/vep_setup.md), ``bias_2015_paths``
writes BIAS's ``required_paths.json`` for the data bundle (BIAS's own script, which also points to
data files shipped inside the BIAS repository; the bundle files differ by annotator),
``bias_2015`` classifies; ACMG codes listed in ``skip_codes`` are written to a ``--skip_list`` file
(PP5/BP6 on/off as an experimental factor).
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from insilico_clinval._runner import Toolbox
from insilico_clinval.adapters import LoadedPipeline, SutFailure
from insilico_clinval.adapters.command import run_black_box, substitute
from insilico_clinval.adapters.parsers import parse_bias_tsv, parse_generic, vep_allele_map
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas.finding import CanonicalFinding
from insilico_clinval.schemas.pipeline import Bias2015Config, CommandTertiary, TertiaryConfig

NIRVANA_ASSEMBLY = "GRCh38"
# BIAS-2015 doc/vep_setup.md (hg38): custom annotation fields BIAS reads from VEP's JSON.
VEP_GNOMAD_FIELDS = "AF%AF_joint%AC%AN"
VEP_CLINVAR_FIELDS = "CLNSIG%CLNREVSTAT%CLNDN%CLNVC%CLNSIGSCV"


def prepare_input(
    toolbox: Toolbox,
    cfg: TertiaryConfig,
    pipeline: LoadedPipeline,
    vcf: Path,
    workdir: Path,
    *,
    sample_id: str,
    manifest: ManifestBuilder | None,
) -> Path:
    """Policy-filtered, biallelic, bgzipped + indexed tertiary input VCF."""
    # absolute: a docker-run bcftools does not share the working directory, so a relative -o
    # (a relative --out) would point nowhere inside the container
    workdir = workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    selected = workdir / "selected.vcf.gz"
    args = ["view", "--no-version", "-O", "z", "-o", str(selected)]
    mounts = [vcf.parent, workdir]
    if cfg.pass_only:
        args += ["-f", "PASS,."]
    if cfg.regions_bed is not None:
        bed = pipeline.resolve(cfg.regions_bed)
        args += ["-T", str(bed)]
        mounts.append(bed.parent)
    toolbox.run(
        "bcftools",
        [*args, str(vcf.resolve())],
        step=f"bcftools view ({sample_id}, tertiary input)",
        manifest=manifest,
        mounts=mounts,
    )
    out = workdir / "tertiary_input.vcf.gz"
    toolbox.run(
        "bcftools",
        ["norm", "--no-version", "-m", "-any", "-O", "z", "-o", str(out), str(selected)],
        step=f"bcftools norm ({sample_id}, tertiary input)",
        manifest=manifest,
        mounts=[workdir],
    )
    toolbox.run(
        "bcftools",
        ["index", "-t", "-f", str(out)],
        step=f"bcftools index ({sample_id}, tertiary input)",
        manifest=manifest,
        mounts=[workdir],
    )
    selected.unlink(missing_ok=True)
    return out


def write_skip_list(codes: Sequence[str], path: Path) -> Path | None:
    if not codes:
        return None
    path.write_text("".join(f"{c}\n" for c in sorted({c.upper() for c in codes})))
    return path


def nirvana_args(
    cfg: Bias2015Config, pipeline: LoadedPipeline, vcf: Path, prefix: Path
) -> list[str]:
    n = cfg.nirvana
    assert n is not None
    return [
        "--cache", str(pipeline.resolve(n.cache)),
        "--sd", str(pipeline.resolve(n.supplementary)),
        "--ref", str(pipeline.resolve(n.reference)),
        "--in", str(vcf.resolve()),
        "--out", str(prefix.resolve()),
    ]  # fmt: skip


def vep_args(cfg: Bias2015Config, pipeline: LoadedPipeline, vcf: Path, out: Path) -> list[str]:
    v = cfg.vep
    assert v is not None
    gnomad = pipeline.resolve(v.gnomad_exomes)
    clinvar = pipeline.resolve(v.clinvar)
    return [
        "--input_file", str(vcf.resolve()),
        "--output_file", str(out.resolve()),
        "--json", "--cache", "--offline",
        "--assembly", NIRVANA_ASSEMBLY,
        "--dir_cache", str(pipeline.resolve(v.cache_dir)),
        "--cache_version", v.cache_version,
        "--fasta", str(pipeline.resolve(v.fasta)),
        "--everything",
        "--fork", str(v.fork),
        "--custom", f"file={gnomad},short_name=gnomAD_exomes,format=vcf,type=exact,coords=0,"
        f"fields={VEP_GNOMAD_FIELDS}",
        "--custom", f"file={clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,"
        f"fields={VEP_CLINVAR_FIELDS}",
        "--plugin", f"AlphaMissense,file={pipeline.resolve(v.alphamissense)}",
        "--plugin", f"REVEL,file={pipeline.resolve(v.revel)}",
        "--force_overwrite", "--no_stats",
    ]  # fmt: skip


def annotate(
    toolbox: Toolbox,
    bias: Bias2015Config,
    pipeline: LoadedPipeline,
    vcf: Path,
    workdir: Path,
    *,
    sample_id: str,
    manifest: ManifestBuilder | None,
) -> Path:
    """Annotator JSON for BIAS-2015 (Nirvana ``<sample>.json.gz`` or VEP ``<sample>.vep.json``)."""
    if bias.annotator == "nirvana":
        assert bias.nirvana is not None
        prefix = workdir / sample_id
        toolbox.run(
            "nirvana",
            nirvana_args(bias, pipeline, vcf, prefix),
            step=f"nirvana ({sample_id})",
            manifest=manifest,
            mounts=[
                workdir,
                pipeline.resolve(bias.nirvana.cache).parent,
                pipeline.resolve(bias.nirvana.supplementary),
                pipeline.resolve(bias.nirvana.reference).parent,
            ],
        )
        annotated = Path(f"{prefix}.json.gz")
    else:
        v = bias.vep
        assert v is not None
        annotated = workdir / f"{sample_id}.vep.json"
        toolbox.run(
            "vep",
            vep_args(bias, pipeline, vcf, annotated),
            step=f"vep ({sample_id})",
            manifest=manifest,
            mounts=[
                workdir,
                vcf.parent,
                pipeline.resolve(v.cache_dir),
                *(
                    pipeline.resolve(p).parent
                    for p in (v.fasta, v.gnomad_exomes, v.clinvar, v.alphamissense, v.revel)
                ),
            ],
        )
        if manifest is not None:
            manifest.add_database("VEP cache", v.cache_version, None)
            manifest.add_database("ClinVar (VEP custom)", v.clinvar_release, None)
    if not annotated.is_file():
        raise SutFailure(f"{bias.annotator} wrote no {annotated.name}")
    return annotated


def run_bias2015(
    toolbox: Toolbox,
    cfg: TertiaryConfig,
    pipeline: LoadedPipeline,
    vcf: Path,
    workdir: Path,
    *,
    sample_id: str,
    manifest: ManifestBuilder | None,
) -> list[CanonicalFinding]:
    bias = cfg.bias2015
    assert bias is not None
    data_dir = pipeline.resolve(bias.data_dir)
    annotated = annotate(
        toolbox, bias, pipeline, vcf, workdir, sample_id=sample_id, manifest=manifest
    )
    paths_json = workdir / "required_paths.json"
    toolbox.run(
        "bias_2015_paths",
        [str(data_dir), "hg38", str(paths_json.resolve()), "--annotator", bias.annotator],
        step=f"bias_2015 required paths ({sample_id})",
        manifest=manifest,
        mounts=[workdir, data_dir],
    )
    out_tsv = workdir / f"{sample_id}.bias_2015.tsv"
    args = [str(annotated.resolve()), str(paths_json.resolve()), str(out_tsv.resolve())]
    args += ["--annotator", bias.annotator]
    skip = write_skip_list(bias.skip_codes, workdir / "skip_list.txt")
    if skip is not None:
        args += ["--skip_list", str(skip.resolve())]
    toolbox.run(
        "bias_2015",
        args,
        step=f"bias_2015 ({sample_id})",
        manifest=manifest,
        mounts=[workdir, data_dir],
    )
    if not out_tsv.is_file():
        raise SutFailure(f"BIAS-2015 wrote no {out_tsv.name}")
    if manifest is not None:
        manifest.add_database("BIAS-2015 data", bias.data_release, None)
    return parse_bias_tsv(
        out_tsv,
        sample_id=sample_id,
        tool="BIAS-2015",
        tool_version=bias.version,
        report_classifications=cfg.report_classifications,
        allele_map=vep_allele_map(vcf) if bias.annotator == "vep" else None,
    )


def run_command_tertiary(
    toolbox: Toolbox,
    cfg: TertiaryConfig,
    vcf: Path,
    workdir: Path,
    *,
    sample_id: str,
    fasta: Path,
    manifest: ManifestBuilder | None,
) -> list[CanonicalFinding]:
    spec: CommandTertiary | None = cfg.command
    assert spec is not None
    values = {
        "sample": sample_id,
        "input_vcf": str(vcf.resolve()),
        "fasta": str(fasta.resolve()),
        "workdir": str(workdir.resolve()),
    }
    output = workdir / substitute(spec.output, values)
    values["output"] = str(output.resolve())
    run_black_box(
        toolbox,
        spec.command,
        values,
        workdir=workdir,
        step=f"{spec.command.tool_name} ({sample_id})",
        manifest=manifest,
        mounts=[vcf.parent, fasta.parent],
    )
    if not output.is_file():
        raise SutFailure(f"{spec.command.tool_name} wrote no {output}")
    return parse_generic(
        output,
        spec.parser,
        sample_id=sample_id,
        tool=spec.command.tool_name,
        tool_version=spec.command.version,
        report_classifications=cfg.report_classifications,
    )
