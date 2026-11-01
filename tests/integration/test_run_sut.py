"""`clinval run-sut` end to end with precomputed / command adapters.

The "pipeline under test" is a stand-in: secondary analysis returns the truth VCF (a perfect
caller), tertiary analysis is tests/fixtures/fake_classifier.py. Needs bcftools, samtools and
Mason (see tests/conftest.py) because run-sut consumes the simulated FASTQ of a real run.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml
from typer.testing import CliRunner

from insilico_clinval.cli import EXIT_SUT, app
from tests.conftest import EREPO_SUITE, FIXTURES, mason_runtime, requires_mason

pytestmark = [pytest.mark.integration, requires_mason]
runner = CliRunner()

COPY = (
    "import shutil, sys, pathlib; out = pathlib.Path(sys.argv[3]); out.mkdir(parents=True, "
    "exist_ok=True); shutil.copy(pathlib.Path(sys.argv[1]) / (sys.argv[2] + '.truth.vcf.gz'), "
    "out / 'calls.vcf.gz')"
)


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    base = tmp_path_factory.mktemp("sut")
    out = base / "run"
    result = runner.invoke(app, ["build-truth", "--suite", str(EREPO_SUITE), "--out", str(out)])
    assert result.exit_code == 0, result.output
    args = ["simulate", "--suite", str(EREPO_SUITE), "--out", str(out)]
    runtime = mason_runtime(base)
    if runtime is not None:
        args += ["--runtime", str(runtime)]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    return out


def _pipeline(path: Path, run: Path, *, secondary: str = "precomputed", fail: bool = False,
              regions: Path | None = None) -> Path:  # fmt: skip
    if secondary == "precomputed":
        sec: dict[str, Any] = {
            "adapter": "precomputed",
            "precomputed": {
                "tool_name": "truth-as-calls",
                "version": "1",
                "samples": {
                    s: {"raw_vcf": str(run / f"{s}.truth.vcf.gz")} for s in ("positive", "negative")
                },
            },
        }
    else:
        sec = {
            "adapter": "command",
            "command": {
                "command": {
                    "tool_name": "copy-caller",
                    "version": "1",
                    "argv": [sys.executable, "-c", COPY, str(run), "{sample}", "{outdir}"],
                },
                "outputs": {"raw_vcf": "calls.vcf.gz"},
            },
        }
    argv = [sys.executable, str(FIXTURES / "fake_classifier.py"), "{input_vcf}", "{output}"]
    if fail:
        argv.append("--fail")
    ter: dict[str, Any] = {
        "adapter": "command",
        "input": "filtered",
        "pass_only": True,
        "command": {
            "command": {"tool_name": "fake-acmg", "version": "0.1", "argv": argv},
            "output": "{sample}.acmg.tsv",
            "parser": {
                "format": "tsv",
                "columns": {
                    "chrom": "CHROM", "pos": "POS", "ref": "REF", "alt": "ALT",
                    "classification": "CLASSIFICATION", "criteria": "CRITERIA", "gene": "GENE",
                },
                "classification_map": {
                    "Pathogenic": "pathogenic", "VUS": "uncertain_significance", "Benign": "benign"
                },
            },
        },
    }  # fmt: skip
    if regions is not None:
        ter["regions_bed"] = str(regions)
    name = f"fake-{secondary}{'-fail' if fail else ''}{'-regions' if regions else ''}"
    path.write_text(
        yaml.safe_dump({"schema_version": 1, "name": name, "secondary": sec, "tertiary": ter})
    )
    return path


def _run_sut(run: Path, pipeline: Path, *extra: str) -> Any:
    return runner.invoke(
        app,
        [
            "run-sut",
            "--suite",
            str(EREPO_SUITE),
            "--out",
            str(run),
            "--pipeline",
            str(pipeline),
            *extra,
        ],
    )


def _findings(run: Path, name: str, sample: str) -> dict[str, Any]:
    return dict(json.loads((run / "sut" / name / f"{sample}.findings.json").read_text()))


def test_run_sut_precomputed_and_command(run_dir: Path, tmp_path: Path) -> None:
    pipeline = _pipeline(tmp_path / "p.yaml", run_dir)
    result = _run_sut(run_dir, pipeline)
    assert result.exit_code == 0, result.output
    pos = _findings(run_dir, "fake-precomputed", "positive")
    assert pos["tool"] == "fake-acmg" and pos["tool_version"] == "0.1"
    planted = json.loads((run_dir / "positive.expected_outcomes.json").read_text())["outcomes"]
    by_pos = {(f["chrom"], f["pos"]): f for f in pos["findings"]}
    for outcome in planted:
        _, chrom, p, _, _ = outcome["variant_id"].split(":")
        finding = by_pos[(chrom, int(p))]
        assert finding["reported"] is True and finding["classification"] == "pathogenic"
        assert finding["acmg_criteria"] == ["PM2", "PVS1"]
    background = [f for f in pos["findings"] if f["classification"] == "uncertain_significance"]
    assert background and not any(f["reported"] for f in background)
    neg = _findings(run_dir, "fake-precomputed", "negative")
    assert {f["classification"] for f in neg["findings"] if f["acmg_criteria"]} == {"benign"}
    secondary = json.loads((run_dir / "sut/fake-precomputed/positive.secondary.json").read_text())
    assert secondary["raw_vcf"] == "positive.truth.vcf.gz" and secondary["adapter"] == "precomputed"
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    assert "pipeline:fake-precomputed" in {i["role"] for i in manifest["inputs"]}
    roles = {o["role"] for o in manifest["outputs"]}
    assert {"secondary:fake-precomputed:positive", "findings:fake-precomputed:negative"} <= roles
    tools = {t["name"]: t["version"] for t in manifest["tools"]}
    assert tools["fake-acmg"] == "0.1" and tools["truth-as-calls"] == "1"
    # tertiary stage alone reuses the stored secondary result
    result = _run_sut(run_dir, pipeline, "--stage", "tertiary", "--sample", "positive")
    assert result.exit_code == 0, result.output


def test_run_sut_command_secondary_and_regions(run_dir: Path, tmp_path: Path) -> None:
    bed = tmp_path / "chr1.bed"
    bed.write_text("chr1\t0\t10000\n")
    pipeline = _pipeline(tmp_path / "c.yaml", run_dir, secondary="command", regions=bed)
    result = _run_sut(run_dir, pipeline, "--sample", "positive")
    assert result.exit_code == 0, result.output
    pos = _findings(run_dir, "fake-command-regions", "positive")
    assert pos["findings"] and {f["chrom"] for f in pos["findings"]} == {"chr1"}
    secondary = json.loads(
        (run_dir / "sut/fake-command-regions/positive.secondary.json").read_text()
    )
    assert secondary["raw_vcf"] == "sut/fake-command-regions/secondary/positive/calls.vcf.gz"


def test_failing_system_under_test(run_dir: Path, tmp_path: Path) -> None:
    pipeline = _pipeline(tmp_path / "f.yaml", run_dir, fail=True)
    result = _run_sut(run_dir, pipeline, "--sample", "positive")
    assert result.exit_code == EXIT_SUT
    assert "simulated crash" in result.output
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    assert any("fake-precomputed-fail failed" in w for w in manifest["warnings"])
