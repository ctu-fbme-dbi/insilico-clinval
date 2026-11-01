"""Fidelity experiment: config, parsers, gaps, planting summary, argv and CLI."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import yaml
from pydantic import ValidationError
from typer.testing import CliRunner

from insilico_clinval._runner import Toolbox
from insilico_clinval.cli import EXIT_INPUT, app
from insilico_clinval.fidelity import load_fidelity, missing_files
from insilico_clinval.fidelity.planting import (
    collect_planting,
    read_planting_tsv,
    summarise_planting,
    write_planting_tsv,
)
from insilico_clinval.fidelity.run import (
    bcftools_stats_args,
    happy_args,
    recall_gaps,
    samtools_stats_args,
)
from insilico_clinval.fidelity.stats import (
    StatsParseError,
    coverage_metrics,
    parse_bcftools_stats,
    parse_happy_extended,
    parse_samtools_stats,
    region_length,
    wilson,
    write_samtools_regions,
)
from insilico_clinval.schemas.fidelity import FidelityConfig
from insilico_clinval.schemas.runtime import RuntimeConfig
from tests.conftest import FIXTURES

runner = CliRunner()

SAMTOOLS_STATS = """\
# This file was produced by samtools stats (1.23) and can be plotted using plot-bamstats
SN\traw total sequences:\t1000\t# excluding supplementary and secondary reads
SN\treads mapped:\t990
SN\treads properly paired:\t980
SN\treads MQ0:\t10\t# mapped and MQ=0
SN\tbases mapped (cigar):\t100000\t# more accurate
SN\terror rate:\t1.500000e-03\t# mismatches / bases mapped (cigar)
SN\tinsert size average:\t400.5
SN\tinsert size standard deviation:\t60.2
IC\t1\t2\t1\t3\t0\t0\t0\t0\t0
IC\t2\t0\t1\t0\t3\t0\t0\t0\t0
COV\t[1-1]\t1\t10
COV\t[10-10]\t10\t20
COV\t[30-30]\t30\t60
"""

BCFTOOLS_STATS = """\
SN\t0\tnumber of samples:\t1
SN\t0\tnumber of SNPs:\t3000
SN\t0\tnumber of indels:\t400
TSTV\t0\t2000\t1000\t2.00\t2000\t1000\t2.00
PSC\t0\tHG002\t0\t1200\t2200\t2000\t1000\t400\t30.0\t0\t0\t0
"""

HAPPY_HEADER = (
    "Type,Subtype,Subset,Filter,Genotype,QQ.Field,QQ,METRIC.Recall,METRIC.Precision,"
    "METRIC.F1_Score,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.TP,QUERY.FP"
)
HAPPY_EXTENDED = "\n".join(
    [
        HAPPY_HEADER,
        "SNP,*,*,ALL,*,QUAL,*,0.99,0.98,0.985,100,99,1,101,99,2",
        "SNP,*,*,PASS,*,QUAL,*,0.95,0.99,0.97,100,95,5,96,95,1",
        "SNP,*,*,PASS,het,QUAL,*,0.95,0.99,0.97,60,57,3,58,57,1",
        "SNP,*,homopolymers,PASS,*,QUAL,*,0.8,0.9,0.85,10,8,2,9,8,1",
        "INDEL,D1_5,*,PASS,*,QUAL,*,0.9,0.9,0.9,10,9,1,10,9,1",
        "INDEL,*,*,PASS,*,QUAL,*,0.9,0.9,0.9,20,18,2,20,18,2",
        "INDEL,*,homopolymers,PASS,*,QUAL,*,nan,nan,nan,0,0,0,0,0,0",
        "SNP,*,*,PASS,*,QUAL,10,0.9,0.9,0.9,100,90,10,91,90,1",
    ]
)


def _config_data(tmp_path: Path) -> dict[str, Any]:
    aln = tmp_path / "sim.bam"
    aln.write_bytes(b"")
    return {
        "schema_version": 1,
        "name": "test",
        "reference": str(FIXTURES / "ref.fa"),
        "confident_bed": str(FIXTURES / "hg002_mini.bed"),
        "stratifications": {"cmrg": str(FIXTURES / "strata_cmrg.bed")},
        "datasets": [
            {
                "id": "simulated",
                "label": "Simulated",
                "alignment": "sim.bam",
                "vcf": str(FIXTURES / "hg002_mini.vcf.gz"),
                "truth_vcf": str(FIXTURES / "hg002_mini.vcf.gz"),
            },
            {
                "id": "real",
                "label": "Real",
                "alignment": "sim.bam",
                "vcf": str(FIXTURES / "hg002_mini.vcf.gz"),
                "truth_vcf": str(FIXTURES / "hg002_mini.vcf.gz"),
            },
        ],
        "baseline": "real",
        "recall_gap": 0.05,
    }


def _write_config(tmp_path: Path, data: dict[str, Any]) -> Path:
    path = tmp_path / "fidelity.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


# -- configuration ----------------------------------------------------------------------


def test_config_validation(tmp_path: Path) -> None:
    data = _config_data(tmp_path)
    FidelityConfig.model_validate(data)
    with pytest.raises(ValidationError, match="Extra inputs"):
        FidelityConfig.model_validate({**data, "unknown": 1})
    with pytest.raises(ValidationError, match="at least 2"):
        FidelityConfig.model_validate({**data, "datasets": data["datasets"][:1]})
    dup = [data["datasets"][0], data["datasets"][0]]
    with pytest.raises(ValidationError, match="unique"):
        FidelityConfig.model_validate({**data, "datasets": dup, "baseline": "simulated"})
    with pytest.raises(ValidationError, match="not a dataset id"):
        FidelityConfig.model_validate({**data, "baseline": "hg002"})


def test_config_paths_resolve_relative_to_file(tmp_path: Path) -> None:
    fid = load_fidelity(_write_config(tmp_path, _config_data(tmp_path)))
    assert fid.resolve("sim.bam") == (tmp_path / "sim.bam").resolve()
    assert fid.qc_regions == fid.confident_bed
    assert missing_files(fid) == []
    data = _config_data(tmp_path)
    data["datasets"][1]["alignment"] = "absent.cram"
    fid = load_fidelity(_write_config(tmp_path, data))
    assert missing_files(fid) == [f"alignment:real: {(tmp_path / 'absent.cram').resolve()}"]


# -- parsers ------------------------------------------------------------------------


def test_wilson() -> None:
    p = wilson(95, 100)
    assert p.value == 0.95
    assert p.lower == pytest.approx(0.8882, abs=1e-4)
    assert p.upper == pytest.approx(0.9785, abs=1e-4)
    assert wilson(0, 0).value is None
    assert wilson(10, 10).upper == 1.0


def test_coverage_counts_uncovered_bases() -> None:
    # 100 bp: 10 at 0x, 10 at 1x, 20 at 10x, 60 at 30x -> 20th percentile = 1x
    mean, ge20, fold80 = coverage_metrics({1: 10, 10: 20, 30: 60}, 100)
    assert mean == pytest.approx((10 + 200 + 1800) / 100)
    assert ge20 == 0.6
    assert fold80 == pytest.approx(20.1)
    with pytest.raises(StatsParseError, match="exceeds"):
        coverage_metrics({1: 200}, 100)


def test_parse_samtools_stats() -> None:
    qc = parse_samtools_stats(SAMTOOLS_STATS, dataset="sim", region_bp=100)
    assert qc.reads_total == 1000
    assert qc.reads_mapped == 0.99
    assert qc.reads_mq0 == 0.01
    assert qc.insert_size_mean == 400.5
    assert qc.error_rate == 0.0015
    assert qc.indel_rate == pytest.approx(10 / 100000)
    assert qc.coverage_ge20 == 0.6
    with pytest.raises(StatsParseError, match="raw total sequences"):
        parse_samtools_stats("SN\treads mapped:\t1\n", dataset="sim", region_bp=100)


def test_parse_bcftools_stats() -> None:
    qc = parse_bcftools_stats(BCFTOOLS_STATS, dataset="real")
    assert (qc.snvs, qc.indels, qc.ti_tv, qc.het_hom) == (3000, 400, 2.0, 1.8333)


def test_parse_happy_extended() -> None:
    rows = parse_happy_extended(HAPPY_EXTENDED, dataset="sim", pass_only=True)
    keys = [(r.type, r.stratum, r.subtype) for r in rows]
    assert keys == [
        ("INDEL", "*", "*"),
        ("INDEL", "*", "D1_5"),
        ("INDEL", "homopolymers", "*"),
        ("SNP", "*", "*"),
        ("SNP", "homopolymers", "*"),
    ]
    snp = rows[3]
    assert (snp.recall.k, snp.recall.n, snp.precision.k, snp.precision.n) == (95, 100, 95, 96)
    assert rows[2].recall.value is None and rows[2].f1 is None
    all_rows = parse_happy_extended(HAPPY_EXTENDED, dataset="sim", pass_only=False)
    assert [(r.recall.k, r.recall.n) for r in all_rows] == [(99, 100)]
    with pytest.raises(StatsParseError, match="missing columns"):
        parse_happy_extended("Type,Subset\nSNP,*\n", dataset="sim", pass_only=True)


def test_samtools_regions_are_one_based_and_merged(tmp_path: Path) -> None:
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t10\t20\nchr1\t15\t30\nchr2\t0\t5\n")
    regions = tmp_path / "r.tsv"
    assert write_samtools_regions(bed, regions) == 25
    assert regions.read_text() == "chr1\t11\t30\nchr2\t1\t5\n"
    assert region_length(regions) == 25


# -- gaps -------------------------------------------------------------------------------


def test_recall_gaps() -> None:
    sim = parse_happy_extended(HAPPY_EXTENDED, dataset="simulated", pass_only=True)
    real_text = HAPPY_EXTENDED.replace(",95,5,96,95,1", ",85,15,86,85,1")
    real = parse_happy_extended(real_text, dataset="real", pass_only=True)
    gaps = recall_gaps(sim + real, "real", 0.05)
    by_key = {(g.type, g.stratum): g for g in gaps}
    assert set(by_key) == {("INDEL", "*"), ("INDEL", "homopolymers"), ("SNP", "*"),
                           ("SNP", "homopolymers")}  # fmt: skip
    assert by_key[("SNP", "*")].gap == pytest.approx(0.10)
    assert by_key[("SNP", "*")].flagged
    assert by_key[("SNP", "homopolymers")].gap == 0 and not by_key[("SNP", "homopolymers")].flagged
    assert by_key[("INDEL", "homopolymers")].gap is None
    assert not any(g.flagged for g in recall_gaps(sim + real, "real", None))


# -- planting ---------------------------------------------------------------------------


def _planting_check(sample: str, statuses: list[tuple[str, str, float | None, int]]) -> str:
    variants = []
    for i, (zyg, status, vaf, alt) in enumerate(statuses):
        variants.append(
            {
                "variant_id": f"curated:chr1:{100 + i}:A:G",
                "chrom": "chr1",
                "pos": 100 + i,
                "ref": "A",
                "alt": "G",
                "zygosity": zyg,
                "carrying_haplotypes": [1],
                "expected_vaf": 0.5 if zyg == "het" else 1.0,
                "spanning_reads": 30,
                "alt_reads": alt,
                "ref_reads": 30 - alt,
                "other_reads": 0,
                "vaf": vaf,
                "by_haplotype": [],
                "status": status,
            }
        )
    data = {
        "schema_version": 1,
        "sample_id": sample,
        "min_alt_reads": 1,
        "window_pad": 10,
        "complete": all(s == "ok" for _, s, _, _ in statuses),
        "variants": variants,
    }
    return json.dumps(data)


def _planting_runs(tmp_path: Path) -> list[Path]:
    a, b = tmp_path / "suite_a", tmp_path / "suite_b"
    a.mkdir()
    b.mkdir()
    (a / "positive.planting_check.json").write_text(
        _planting_check("positive", [("het", "ok", 0.5, 15), ("hemi", "ok", 1.0, 14)])
    )
    (a / "negative.planting_check.json").write_text(
        _planting_check("negative", [("het", "ok", 0.4, 12)])
    )
    (b / "positive.planting_check.json").write_text(
        _planting_check("positive", [("het", "no_reads", 0.0, 0)])
    )
    return [a, b]


def test_planting_summary(tmp_path: Path) -> None:
    rows = collect_planting(_planting_runs(tmp_path))
    assert [(r.run, r.sample_id) for r in rows] == [
        ("suite_a", "negative"),
        ("suite_a", "positive"),
        ("suite_a", "positive"),
        ("suite_b", "positive"),
    ]
    tsv = write_planting_tsv(rows, tmp_path / "planting.tsv")
    assert read_planting_tsv(tsv) == rows
    s = summarise_planting(rows)
    assert (s.samples, s.variants, s.complete.k, s.complete.n) == (3, 4, 3, 4)
    assert s.by_status == {"ok": 3, "no_reads": 1}
    het = next(z for z in s.by_zygosity if z.zygosity.value == "het")
    assert (het.variants, het.vaf_median, het.vaf_min, het.alt_reads_min) == (3, 0.4, 0.0, 0)


# -- argv through the Toolbox -----------------------------------------------------------


def test_argv_and_docker_mounts(tmp_path: Path) -> None:
    data = _config_data(tmp_path)
    data["happy"] = {"vcfeval_template": str(tmp_path / "ref.sdf")}
    fid = load_fidelity(_write_config(tmp_path, data))
    d = fid.dataset("simulated")
    regions = tmp_path / "qc" / "regions.tsv"
    args = samtools_stats_args(fid, d, regions)
    assert args[:3] == ["stats", "--reference", str(FIXTURES / "ref.fa")]
    assert args[-3:] == ["-t", str(regions), str((tmp_path / "sim.bam").resolve())]
    bargs = bcftools_stats_args(fid, d)
    assert ["-T", str(FIXTURES / "hg002_mini.bed")] == bargs[5:7]
    assert bargs[-3:-1] == ["-f", "PASS"]
    prefix = tmp_path / "happy" / "simulated" / "simulated"
    hargs = happy_args(fid, d, prefix, tmp_path / "strata.tsv")
    for flag in ("--engine-vcfeval-template", "--pass-only", "--stratification", "-f", "-r"):
        assert flag in hargs
    assert "--gender" not in hargs
    toolbox = Toolbox(
        RuntimeConfig.model_validate(
            {"tools": {"hap.py": {"runner": "docker", "image": "h@sha256:1"}}}
        )
    )
    argv = toolbox.command("hap.py", hargs, mounts=[FIXTURES, tmp_path], cwd=prefix.parent)
    mounts = [argv[i + 1] for i, a in enumerate(argv) if a == "-v"]
    assert mounts == sorted(
        [f"{FIXTURES}:{FIXTURES}", f"{tmp_path.resolve()}:{tmp_path.resolve()}"]
    )
    assert argv[argv.index("h@sha256:1") + 1] == "hap.py"


# -- CLI --------------------------------------------------------------------------------


def test_cli_planting_step_and_assembly_from_raw_outputs(tmp_path: Path) -> None:
    data = _config_data(tmp_path)
    data["planting"] = [str(p) for p in _planting_runs(tmp_path)]
    config = _write_config(tmp_path, data)
    out = tmp_path / "out"
    (out / "qc").mkdir(parents=True)
    (out / "qc" / "regions.samtools.tsv").write_text("chr1\t1\t100\n")
    (out / "qc" / "simulated.samtools_stats.txt").write_text(SAMTOOLS_STATS)
    (out / "qc" / "simulated.bcftools_stats.txt").write_text(BCFTOOLS_STATS)
    for ds in ("simulated", "real"):
        (out / "happy" / ds).mkdir(parents=True)
        (out / "happy" / ds / f"{ds}.extended.csv").write_text(HAPPY_EXTENDED)

    result = runner.invoke(
        app, ["fidelity", "--config", str(config), "--out", str(out), "--step", "planting"]
    )
    assert result.exit_code == 0, result.output
    res = json.loads((out / "fidelity.json").read_text())
    assert res["planting"]["variants"] == 4
    assert [q["dataset"] for q in res["read_qc"]] == ["simulated"]
    assert res["missing"] == ["qc/real.samtools_stats.txt", "qc/real.bcftools_stats.txt"]
    assert {g["gap"] for g in res["gaps"]} == {0.0, None}
    assert (out / "fidelity_calling.tsv").read_text().startswith("dataset\ttype\tstratum")
    manifest = json.loads((out / "run_manifest.json").read_text())
    roles = {i["role"] for i in manifest["inputs"]}
    assert roles == {"planting_check:suite_a", "planting_check:suite_b"}
    assert "fidelity: raw output not found: qc/real.samtools_stats.txt" in manifest["warnings"]
    first = (out / "fidelity.json").read_bytes()
    result = runner.invoke(
        app, ["fidelity", "--config", str(config), "--out", str(out), "--step", "planting"]
    )
    assert result.exit_code == 0, result.output
    assert (out / "fidelity.json").read_bytes() == first


def test_cli_errors(tmp_path: Path) -> None:
    data = _config_data(tmp_path)
    config = _write_config(tmp_path, data)
    out = tmp_path / "out"
    result = runner.invoke(
        app, ["fidelity", "--config", str(config), "--out", str(out), "--step", "align"]
    )
    assert result.exit_code == EXIT_INPUT
    assert "unknown --step" in result.output
    data["datasets"][0]["vcf"] = "absent.vcf.gz"
    result = runner.invoke(
        app, ["fidelity", "--config", str(_write_config(tmp_path, data)), "--out", str(out)]
    )
    assert result.exit_code == EXIT_INPUT
    assert "missing inputs" in result.output
    data = _config_data(tmp_path)
    data["planting"] = [str(tmp_path)]
    result = runner.invoke(
        app,
        ["fidelity", "--config", str(_write_config(tmp_path, data)), "--out", str(out),
         "--step", "planting"],
    )  # fmt: skip
    assert result.exit_code == EXIT_INPUT
    assert "no *.planting_check.json" in result.output


def test_happy_gender(tmp_path: Path) -> None:
    """hap.py gets --gender only when configured."""
    data = _config_data(tmp_path)
    data["happy"] = {"gender": "male"}
    fid = load_fidelity(_write_config(tmp_path, data))
    d = fid.dataset("simulated")
    hargs = happy_args(fid, d, tmp_path / "p", None)
    assert hargs[-2:] == ["--gender", "male"]
