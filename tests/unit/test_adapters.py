"""Adapters: pipeline config, sarek, black-box commands, parsers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml
from pydantic import ValidationError

from insilico_clinval._exec import ToolError
from insilico_clinval._runner import Toolbox, version_line
from insilico_clinval.adapters import (
    PipelineConfigError,
    SampleInput,
    SutFailure,
    load_pipeline,
    missing_files,
    referenced_files,
)
from insilico_clinval.adapters import sarek as sarek_adapter
from insilico_clinval.adapters.command import substitute
from insilico_clinval.adapters.parsers import (
    bias_criteria,
    parse_bias_tsv,
    parse_generic,
    vep_allele_map,
)
from insilico_clinval.adapters.tertiary import (
    nirvana_args,
    prepare_input,
    vep_args,
    write_skip_list,
)
from insilico_clinval.schemas.enums import Classification, Sex
from insilico_clinval.schemas.pipeline import ParserSpec, PipelineConfig, TertiaryConfig

REPO = Path(__file__).resolve().parents[2]
P_LP = [Classification.PATHOGENIC, Classification.LIKELY_PATHOGENIC]

SAREK: dict[str, Any] = {
    "schema_version": 1,
    "name": "p1",
    "secondary": {
        "adapter": "sarek",
        "sarek": {
            "release": "3.10.0",
            "tools": ["haplotypecaller"],
            "skip_tools": ["baserecalibrator"],
            "reference": {"fasta": "ref/ref.fa", "fasta_fai": "ref/ref.fa.fai"},
            "params": {"bwa_mem_args": "-K 100000000"},
            "resource_limits": {"cpus": 4, "memory": "6.GB"},
        },
    },
    "tertiary": {"adapter": "none"},
}


def _write(tmp_path: Path, data: dict[str, Any], name: str = "pipeline.yaml") -> Path:
    path = tmp_path / name
    path.write_text(yaml.safe_dump(data))
    return path


# -- config validation -----------------------------------------------------------------


def test_pipeline_config_validation() -> None:
    PipelineConfig.model_validate(SAREK)
    bad = json.loads(json.dumps(SAREK))
    bad["secondary"]["sarek"] = None
    with pytest.raises(ValidationError, match="needs a `sarek` section"):
        PipelineConfig.model_validate(bad)
    bad = json.loads(json.dumps(SAREK))
    bad["tertiary"]["bias2015"] = {"version": "v3.0.0"}
    with pytest.raises(ValidationError):
        PipelineConfig.model_validate(bad)
    bad = json.loads(json.dumps(SAREK))
    bad["secondary"]["sarek"]["reference"] = {"genome": "GATK.GRCh38", "fasta": "x.fa"}
    with pytest.raises(ValidationError, match="exactly one of genome"):
        PipelineConfig.model_validate(bad)
    bad = json.loads(json.dumps(SAREK))
    bad["secondary"]["sarek"]["typo"] = 1
    with pytest.raises(ValidationError, match="Extra inputs"):
        PipelineConfig.model_validate(bad)
    with pytest.raises(ValidationError, match=r"parser\.columns misses"):
        ParserSpec.model_validate(
            {"format": "tsv", "columns": {"chrom": "C"}, "classification_map": {}}
        )


def test_load_pipeline_and_missing_files(tmp_path: Path) -> None:
    path = _write(tmp_path, SAREK)
    loaded = load_pipeline(path)
    assert loaded.resolve("ref/ref.fa") == (tmp_path / "ref/ref.fa").resolve()
    assert [m.split(":")[0] for m in missing_files(loaded)] == ["sarek_fasta", "sarek_fasta_fai"]
    (tmp_path / "bad.yaml").write_text("schema_version: 1\nname: x\n")
    with pytest.raises(PipelineConfigError, match="invalid pipeline config"):
        load_pipeline(tmp_path / "bad.yaml")


def test_repository_pipeline_configs_are_valid() -> None:
    for path in sorted((REPO / "configs/adapters").glob("*.yaml")):
        load_pipeline(path)
    load_pipeline(REPO / "tests/fixtures/pipeline_sarek_fixture.yaml")


# -- sarek -----------------------------------------------------------------------------


def _samples(tmp_path: Path) -> list[SampleInput]:
    return [
        SampleInput("positive", Sex.MALE, tmp_path / "p.R1.fq.gz", tmp_path / "p.R2.fq.gz"),
        SampleInput("negative", Sex.FEMALE, tmp_path / "n.R1.fq.gz", tmp_path / "n.R2.fq.gz"),
    ]


def test_sarek_prepare(tmp_path: Path) -> None:
    pipeline = load_pipeline(_write(tmp_path, SAREK))
    assert pipeline.config.secondary.sarek is not None
    work = tmp_path / "secondary"
    args = sarek_adapter.prepare(
        pipeline.config.secondary.sarek, pipeline, _samples(tmp_path), work
    )
    sheet = (work / "samplesheet.csv").read_text().splitlines()
    assert sheet[0] == "patient,sex,status,sample,lane,fastq_1,fastq_2"
    assert sheet[1].startswith("positive,XY,0,positive,lane_1,")
    assert sheet[2].startswith("negative,XX,0,negative,lane_1,")
    params = yaml.safe_load((work / "params.yaml").read_text())
    assert params["igenomes_ignore"] is True and "genome" not in params
    assert params["fasta"] == str((tmp_path / "ref/ref.fa").resolve())
    assert params["tools"] == "haplotypecaller" and params["skip_tools"] == "baserecalibrator"
    assert params["bwa_mem_args"] == "-K 100000000"
    assert params["outdir"] == str((work / "results").resolve())
    config = (work / "resources.config").read_text()
    assert "resourceLimits = [ cpus: 4, memory: '6.GB' ]" in config
    assert args[:6] == ["run", "nf-core/sarek", "-r", "3.10.0", "-profile", "docker"]
    assert "-params-file" in args and "-c" in args and "-work-dir" in args


def test_sarek_reserved_params(tmp_path: Path) -> None:
    data = json.loads(json.dumps(SAREK))
    data["secondary"]["sarek"]["params"] = {"outdir": "/x"}
    pipeline = load_pipeline(_write(tmp_path, data))
    assert pipeline.config.secondary.sarek is not None
    with pytest.raises(PipelineConfigError, match="must not set"):
        sarek_adapter.prepare(
            pipeline.config.secondary.sarek, pipeline, _samples(tmp_path), tmp_path
        )


def test_sarek_collect_outputs(tmp_path: Path) -> None:
    results = tmp_path / "results"
    calls = results / "variant_calling/haplotypecaller/s1"
    with pytest.raises(SutFailure, match="no haplotypecaller VCF"):
        sarek_adapter.collect_outputs(results, "s1")
    calls.mkdir(parents=True)
    (calls / "s1.haplotypecaller.vcf.gz").write_text("x")
    md = results / "preprocessing/markduplicates/s1"
    md.mkdir(parents=True)
    (md / "s1.md.cram").write_text("x")
    out = sarek_adapter.collect_outputs(results, "s1")
    assert out.filtered_vcf is None and out.alignment == md / "s1.md.cram"
    recal = results / "preprocessing/recalibrated/s1"
    recal.mkdir(parents=True)
    (recal / "s1.recal.cram").write_text("x")
    (calls / "s1.haplotypecaller.filtered.vcf.gz").write_text("x")
    out = sarek_adapter.collect_outputs(results, "s1")
    assert out.alignment == recal / "s1.recal.cram"
    assert out.filtered_vcf == calls / "s1.haplotypecaller.filtered.vcf.gz"


def test_sarek_adopt_external_run(tmp_path: Path) -> None:
    data = json.loads(json.dumps(SAREK))
    data["secondary"]["sarek"].update(
        {
            "tools": ["deepvariant", "vep"],
            "variant_caller": "deepvariant",
            "reference": {"genome": "GATK.GRCh38"},
            "params": {"aligner": "bwa-mem2", "trim_fastq": True},
        }
    )
    pipeline = load_pipeline(_write(tmp_path, data))
    cfg = pipeline.config.secondary.sarek
    assert cfg is not None
    samples = _samples(tmp_path)
    launch = tmp_path / "launch"
    results = launch / "results"
    with pytest.raises(SutFailure, match=r"params\.yaml not found"):
        sarek_adapter.adopt(cfg, pipeline, samples, results)
    launch.mkdir()
    sarek_adapter.write_samplesheet(samples, launch / "samplesheet.csv")
    used = {
        "input": "x", "outdir": "y", "genome": "GATK.GRCh38", "tools": "deepvariant,vep",
        "skip_tools": "baserecalibrator", "aligner": "bwa-mem2", "trim_fastq": True,
        "igenomes_base": "/site/igenomes",
    }  # fmt: skip
    (launch / "params.yaml").write_text(yaml.safe_dump(used))
    (results / "pipeline_info").mkdir(parents=True)
    versions = results / sarek_adapter.VERSIONS_FILE
    versions.write_text("Workflow:\n    nf-core/sarek: v3.10.0\n")
    evidence = sarek_adapter.adopt(cfg, pipeline, samples, results)
    assert [p.name for p in evidence] == ["params.yaml", "samplesheet.csv", versions.name]
    (launch / "params.yaml").write_text(yaml.safe_dump({**used, "trim_fastq": False}))
    versions.write_text("nf-core/sarek: v3.9.0\n")
    with pytest.raises(SutFailure) as exc:
        sarek_adapter.adopt(cfg, pipeline, samples, results)
    assert "trim_fastq" in str(exc.value) and "release" in str(exc.value)
    dv = results / "variant_calling/deepvariant/positive"
    dv.mkdir(parents=True)
    (dv / "positive.deepvariant.vcf.gz").write_text("x")
    out = sarek_adapter.collect_outputs(results, "positive", "deepvariant")
    assert out.raw_vcf == dv / "positive.deepvariant.vcf.gz" and out.filtered_vcf is None


def test_sarek_deepvariant_haploid(tmp_path: Path) -> None:
    data = json.loads(json.dumps(SAREK))
    data["secondary"]["sarek"]["deepvariant_haploid"] = {"par_bed": "par.bed"}
    with pytest.raises(PipelineConfigError, match="needs `deepvariant`"):
        load_pipeline(_write(tmp_path, data))
    data["secondary"]["sarek"].update({"tools": ["deepvariant"], "variant_caller": "deepvariant"})
    (tmp_path / "par.bed").write_text("chrX\t10000\t2781479\n")
    pipeline = load_pipeline(_write(tmp_path, data))
    cfg = pipeline.config.secondary.sarek
    assert cfg is not None
    text = sarek_adapter.deepvariant_config(cfg, pipeline)
    par = (tmp_path / "par.bed").resolve()
    assert text is not None
    assert f"meta.sex == 'XY' ? \"--haploid_contigs=chrX,chrY --par_regions_bed={par}\"" in text
    assert f'containerOptions = "-v {par.parent}:{par.parent}"' in text
    assert '"--sample_name ${meta.sample}"' in text
    work = tmp_path / "secondary"
    args = sarek_adapter.prepare(cfg, pipeline, _samples(tmp_path), work)
    assert (work / "deepvariant.config").read_text() == text
    assert args.count("-c") == 2 and str((work / "deepvariant.config").resolve()) in args
    # adoption: the launch directory must hold the same deepvariant.config
    results = work / "results"
    (results / "pipeline_info").mkdir(parents=True)
    (results / sarek_adapter.VERSIONS_FILE).write_text("nf-core/sarek: v3.10.0\n")
    evidence = sarek_adapter.adopt(cfg, pipeline, _samples(tmp_path), results)
    assert evidence[-1].name == "deepvariant.config"
    (work / "deepvariant.config").write_text("// other\n")
    with pytest.raises(SutFailure, match=r"deepvariant\.config: missing or different"):
        sarek_adapter.adopt(cfg, pipeline, _samples(tmp_path), results)


# -- BIAS-2015 ---------------------------------------------------------------------------

BIAS_HEADER = (
    "chromosome\tposition\trefAllele\taltAllele\tvariantType\tconsequence\tacmgClassification\t"
    "alleleFreq\thgvsg\thgvsc\thgvsp\taaChange\tgeneName\tpubmedIds\tassociatedDiseases\t"
    "dbSnpids\ttranscript\trationale\n"
)


def _bias_row(chrom: str, pos: int, ref: str, alt: str, label: str, met: dict[str, int]) -> str:
    rationale = {
        "pvs": {"pvs1": [met.get("pvs1", 0), "x"]},
        "pm": {"pm2": [met.get("pm2", 0), ""]},
        "ba": {"ba1": [met.get("ba1", 0), ""]},
    }
    return (
        f"{chrom}\t{pos}\t{ref}\t{alt}\tSNV\tstop_gained\t{label}\t\t\tNM_1:c.1A>T\t\t\tGENEA\t"
        f"\t\t\tNM_1\t{json.dumps(rationale)}\n"
    )


def test_parse_bias_tsv(tmp_path: Path) -> None:
    path = tmp_path / "b.tsv"
    path.write_text(
        BIAS_HEADER
        + _bias_row("chr1", 100, "a", "t", "pathogenic", {"pvs1": 3, "pm2": 1})
        + _bias_row("chr1", 200, "C", "G", "likely pathogenic", {"pm2": 1})
        + _bias_row("chr1", 300, "C", "G", "uncertain", {})
        + _bias_row("chr1", 400, "C", "G", "benign", {"ba1": 1})
    )
    findings = parse_bias_tsv(
        path, sample_id="s", tool="BIAS-2015", tool_version="v3.0.0", report_classifications=P_LP
    )
    assert [(f.pos, f.classification, f.reported) for f in findings] == [
        (100, Classification.PATHOGENIC, True),
        (200, Classification.LIKELY_PATHOGENIC, True),
        (300, Classification.UNCERTAIN_SIGNIFICANCE, False),
        (400, Classification.BENIGN, False),
    ]
    assert findings[0].acmg_criteria == ["PM2", "PVS1"] and findings[0].ref == "A"
    assert findings[0].gene == "GENEA" and findings[0].hgvs_c == "NM_1:c.1A>T"
    assert findings[0].hgvs_p is None
    path.write_text(BIAS_HEADER + _bias_row("chr1", 1, "A", "T", "conflicting", {}))
    with pytest.raises(SutFailure, match="unknown BIAS classification"):
        parse_bias_tsv(path, sample_id="s", tool="t", tool_version="v", report_classifications=P_LP)
    assert bias_criteria("") == []


def test_parse_bias_tsv_vep_alleles(tmp_path: Path) -> None:
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr10,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr10\t100\t.\tCAA\tC\t.\tPASS\t.\n"
        "chr10\t200\t.\tG\tGTA\t.\tPASS\t.\n"
        "chr10\t300\t.\tC\tT\t.\tPASS\t.\n"
        "chr10\t400\t.\tCAAAACAAAAACA\tCAAAACA\t.\tPASS\t.\n"
        "chr10\t500\t.\tCATTCATTTATTTATTTATTT\tCATTTATTTATTT\t.\tPASS\t.\n"
        "chr10\t600\t.\tAAAGAAGAAGAAG\tGAAGAAGAAGAAGAAGAAG\t.\tPASS\t.\n"
    )
    amap = vep_allele_map(vcf)
    assert amap[("chr10", 407, "AAAACA", "-")] == ("chr10", 400, "CAAAACAAAAACA", "CAAAACA")
    assert amap[("chr10", 300, "C", "T")] == ("chr10", 300, "C", "T")
    assert amap[("chr10", 504, "CATTTATT", "-")][1:] == (
        500,
        "CATTCATTTATTTATTTATTT",
        "CATTTATTTATTT",
    )
    assert amap[("chr10", 101, "AA", "-")] == ("chr10", 100, "CAA", "C")
    assert amap[("chr10", 201, "-", "TA")] == ("chr10", 200, "G", "GTA")
    # no shared first base: VEP trims the shared suffix only (chr7:127015001, 2026-10-30)
    assert amap[("chr10", 600, "A", "GAAGAAG")][2:] == ("AAAGAAGAAGAAG", "GAAGAAGAAGAAGAAGAAG")
    path = tmp_path / "bias.tsv"
    path.write_text(
        BIAS_HEADER
        + _bias_row("chr10", 101, "AA", "-", "pathogenic", {"pvs1": 3})
        + _bias_row("chr10", 201, "-", "TA", "benign", {})
        + _bias_row("chr10", 300, "C", "T", "benign", {})
    )
    args: dict[str, Any] = {"sample_id": "s", "tool": "t", "tool_version": "v"}
    found = parse_bias_tsv(path, report_classifications=P_LP, allele_map=amap, **args)
    assert [(f.pos, f.ref, f.alt, f.reported) for f in found] == [
        (100, "CAA", "C", True),
        (200, "G", "GTA", False),
        (300, "C", "T", False),
    ]
    path.write_text(BIAS_HEADER + _bias_row("chr10", 999, "A", "-", "benign", {}))
    with pytest.raises(SutFailure, match="matches no record"):
        parse_bias_tsv(path, report_classifications=P_LP, allele_map=amap, **args)


def test_skip_list_and_nirvana_args(tmp_path: Path) -> None:
    assert write_skip_list([], tmp_path / "skip.txt") is None
    path = write_skip_list(["bp6", "PP5", "PP5"], tmp_path / "skip.txt")
    assert path is not None and path.read_text() == "BP6\nPP5\n"
    data = json.loads(json.dumps(SAREK))
    data["tertiary"] = {
        "adapter": "bias2015",
        "bias2015": {
            "version": "v3.0.0",
            "data_dir": "bias",
            "data_release": "2026.03.01",
            "nirvana": {"cache": "n/Cache/Both", "supplementary": "n/SA", "reference": "n/R.dat"},
            "skip_codes": ["PP5", "BP6"],
        },
    }
    pipeline = load_pipeline(_write(tmp_path, data))
    assert pipeline.config.tertiary.bias2015 is not None
    args = nirvana_args(
        pipeline.config.tertiary.bias2015, pipeline, tmp_path / "in.vcf.gz", tmp_path / "s"
    )
    assert args[0:2] == ["--cache", str((tmp_path / "n/Cache/Both").resolve())]
    assert args[-2:] == ["--out", str((tmp_path / "s").resolve())]


VEP_DATA = {
    "cache_dir": "vep/cache",
    "cache_version": "115",
    "fasta": "vep/cache/genome.fa.gz",
    "gnomad_exomes": "vep/gnomad.vcf.bgz",
    "clinvar": "vep/clinvar.vcf.gz",
    "clinvar_release": "2026-02-01",
    "alphamissense": "vep/am.tsv.gz",
    "revel": "vep/revel.tsv.gz",
}


def test_vep_annotator(tmp_path: Path) -> None:
    data = json.loads(json.dumps(SAREK))
    bias = {"version": "v3.0.0", "data_dir": "bias", "data_release": "2026.03.01"}
    data["tertiary"] = {"adapter": "bias2015", "bias2015": {**bias, "annotator": "vep"}}
    with pytest.raises(ValidationError, match="adapter vep needs a `vep` section"):
        PipelineConfig.model_validate(data)
    nirvana = {"cache": "n/C", "supplementary": "n/SA", "reference": "n/R.dat"}
    data["tertiary"]["bias2015"] = {**bias, "annotator": "vep", "vep": VEP_DATA, "nirvana": nirvana}
    with pytest.raises(ValidationError, match="section `nirvana` given but adapter is vep"):
        PipelineConfig.model_validate(data)
    data["tertiary"]["bias2015"] = {**bias, "annotator": "vep", "vep": VEP_DATA}
    pipeline = load_pipeline(_write(tmp_path, data))
    b = pipeline.config.tertiary.bias2015
    assert b is not None
    roles = [role for role, _ in referenced_files(pipeline)]
    assert "vep_gnomad_exomes" in roles and "nirvana_cache" not in roles
    args = vep_args(b, pipeline, tmp_path / "in.vcf.gz", tmp_path / "s.vep.json")
    assert args[args.index("--output_file") + 1] == str((tmp_path / "s.vep.json").resolve())
    assert (
        "--json" in args
        and "--offline" in args
        and args[args.index("--cache_version") + 1] == "115"
    )
    customs = [args[i + 1] for i, a in enumerate(args) if a == "--custom"]
    gnomad = str((tmp_path / "vep/gnomad.vcf.bgz").resolve())
    assert customs[0].startswith(f"file={gnomad},short_name=gnomAD_exomes,")
    assert customs[0].endswith("fields=AF%AF_joint%AC%AN")
    assert "short_name=ClinVar" in customs[1] and "CLNREVSTAT" in customs[1]
    plugins = [args[i + 1] for i, a in enumerate(args) if a == "--plugin"]
    assert [p.split(",")[0] for p in plugins] == ["AlphaMissense", "REVEL"]


# -- generic black box -------------------------------------------------------------------


def test_substitute() -> None:
    assert substitute("{sample}.tsv", {"sample": "s1"}) == "s1.tsv"
    with pytest.raises(PipelineConfigError, match="unknown placeholder"):
        substitute("{vcf}", {"input_vcf": "x"})


SPEC: dict[str, Any] = {
    "format": "tsv",
    "columns": {
        "chrom": "CHROM", "pos": "POS", "ref": "REF", "alt": "ALT",
        "classification": "CLASS", "criteria": "CODES",
    },
    "classification_map": {"P": "pathogenic", "VUS": "uncertain_significance", "B": "benign"},
}  # fmt: skip


def test_parse_generic_tsv_and_json(tmp_path: Path) -> None:
    path = tmp_path / "o.tsv"
    path.write_text(
        "CHROM\tPOS\tREF\tALT\tCLASS\tCODES\nchr1\t5\tA\tT\tP\tpvs1, PM2\nchr1\t9\tA\tT\t\t\n"
    )
    spec = ParserSpec.model_validate(SPEC)
    findings = parse_generic(
        path, spec, sample_id="s", tool="x", tool_version="1", report_classifications=P_LP
    )
    assert [(f.classification, f.reported, f.acmg_criteria) for f in findings] == [
        (Classification.PATHOGENIC, True, ["PM2", "PVS1"]),
        (None, False, []),
    ]
    js = tmp_path / "o.json"
    js.write_text(
        json.dumps(
            {"variants": [{"c": "chr2", "p": 7, "r": "G", "a": "C", "k": "B", "rep": "yes"}]}
        )
    )
    spec_json = ParserSpec.model_validate(
        {
            "format": "json",
            "records": "variants",
            "columns": {
                "chrom": "c",
                "pos": "p",
                "ref": "r",
                "alt": "a",
                "classification": "k",
                "reported": "rep",
            },
            "classification_map": {"B": "benign"},
        }
    )
    [f] = parse_generic(
        js, spec_json, sample_id="s", tool="x", tool_version="1", report_classifications=P_LP
    )
    assert (f.chrom, f.pos, f.classification, f.reported) == (
        "chr2",
        7,
        Classification.BENIGN,
        True,
    )


def test_parse_generic_errors(tmp_path: Path) -> None:
    spec = ParserSpec.model_validate(SPEC)
    path = tmp_path / "o.tsv"
    path.write_text("CHROM\tPOS\tREF\tALT\tCLASS\tCODES\nchr1\t5\tA\tT\tLP?\t\n")
    with pytest.raises(SutFailure, match="not in classification_map"):
        parse_generic(
            path, spec, sample_id="s", tool="x", tool_version="1", report_classifications=P_LP
        )
    path.write_text("CHROM\tPOS\tREF\n")
    with pytest.raises(SutFailure, match="missing columns"):
        parse_generic(
            path, spec, sample_id="s", tool="x", tool_version="1", report_classifications=P_LP
        )


def test_black_box_commands(tmp_path: Path) -> None:
    box = Toolbox()
    box.register_black_box("closed-acmg", runner="docker", image="acmg@sha256:ab", version="2.1")
    argv = box.command(
        "closed-acmg", ["classify", "--in", "x.vcf"], mounts=[tmp_path], cwd=tmp_path
    )
    assert argv[-4:] == ["acmg@sha256:ab", "classify", "--in", "x.vcf"]
    assert box.version("closed-acmg") == "2.1"
    box.register_black_box("local-tool", runner="local", image=None, version="1")
    assert box.command("local-tool", ["a", "b"]) == ["a", "b"]
    with pytest.raises(ToolError, match="clashes"):
        box.register_black_box("bcftools", runner="local", image=None, version="1")


def test_version_line() -> None:
    banner = "\n      N E X T F L O W\n      version 26.04.6 build 0\n"
    assert version_line(banner, banner="version") == "version 26.04.6 build 0"
    vep_help = "# ENSEMBL VARIANT EFFECT PREDICTOR #\nVersions:\n  ensembl-vep : 115.0\n"
    assert version_line(vep_help, banner="ensembl-vep") == "ensembl-vep : 115.0"
    assert version_line("bcftools 1.23\nUsing htslib\nGNU GPL version 3") == "bcftools 1.23"
    assert version_line("") == "unknown"


class _RecordingToolbox(Toolbox):
    def __init__(self) -> None:
        super().__init__()
        self.calls: list[tuple[str, list[str], list[Path]]] = []

    def run(self, tool: str, args: Any, **kwargs: Any) -> Any:  # type: ignore[override]
        self.calls.append((tool, list(args), list(kwargs.get("mounts", []))))


def test_prepare_input_paths_absolute(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # a relative --out must not reach a docker-run bcftools as a relative path (2026-10-30)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "in").mkdir()
    box = _RecordingToolbox()
    out = prepare_input(
        box,
        TertiaryConfig(adapter="none"),
        None,  # type: ignore[arg-type]  # only read for regions_bed
        tmp_path / "in" / "s.vcf.gz",
        Path("runs/x/tertiary/s"),
        sample_id="s",
        manifest=None,
    )
    assert out.is_absolute()
    assert len(box.calls) == 3
    for _, args, mounts in box.calls:
        assert all(Path(a).is_absolute() for a in args if a.endswith(".vcf.gz"))
        assert all(m.is_absolute() for m in mounts)
