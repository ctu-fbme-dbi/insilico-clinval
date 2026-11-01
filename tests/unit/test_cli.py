"""CLI exit codes and outputs without external tools."""

from __future__ import annotations

import json
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from insilico_clinval.cli import (
    EXIT_INPUT,
    EXIT_SELECTION,
    EXIT_TRUTH,
    LOG_NAME,
    PROGRESS_LOGGERS,
    _run_logging,
    app,
)
from insilico_clinval.provenance import sha256_file
from tests.conftest import EREPO_SUITE, FIXTURES, SUITE, rule

runner = CliRunner()


def _run_dir(suite: Path) -> Path:
    return suite.parent / f"run_{suite.stem}"


def test_validate_suite_ok() -> None:
    result = runner.invoke(app, ["validate-suite", str(SUITE)])
    assert result.exit_code == 0, result.output
    assert "6 rules" in result.output


def test_validate_suite_missing_input(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    suite_data["sources"]["clinvar"]["vcf"] = "/nonexistent/clinvar.vcf.gz"
    result = runner.invoke(app, ["validate-suite", str(write_suite(suite_data))])
    assert result.exit_code == EXIT_INPUT
    assert "missing required inputs" in result.output


def test_select_writes_outputs_and_manifest(tmp_path: Path) -> None:
    out = tmp_path / "run"
    result = runner.invoke(app, ["select", "--suite", str(SUITE), "--out", str(out)])
    assert result.exit_code == 0, result.output
    for name in ("variants.selected.json", "variants.selected.tsv", "selection_report.md"):
        assert (out / name).is_file()
    manifest = json.loads((out / "run_manifest.json").read_text())
    inputs = {i["role"]: i for i in manifest["inputs"]}
    assert inputs["clinvar_vcf"]["sha256"] == sha256_file(FIXTURES / "clinvar_mini.vcf")
    assert inputs["reference_fasta"]["sha256"] == sha256_file(FIXTURES / "ref.fa")
    assert {d["name"] for d in manifest["databases"]} == {"ClinVar", "GIAB"}
    assert manifest["seed"] == 20260928
    assert [s["name"] for s in manifest["steps"]] == ["select"]


def test_unsatisfied_required_rule_fails(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    """A required rule that cannot be satisfied -> non-zero exit, explained in the report."""
    rule(suite_data, "indel_ladder")["required"] = True
    suite = write_suite(suite_data)
    result = runner.invoke(app, ["select", "--suite", str(suite)])
    assert result.exit_code == EXIT_SELECTION
    assert "indel_ladder" in result.output
    report = (_run_dir(suite) / "selection_report.md").read_text()
    assert "## FAILED" in report
    assert "required rule indel_ladder is unsatisfied" in report
    assert "21-50bp: requested 2, selected 1" in report


def test_empty_sample_fails(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    """A sample without variants -> non-zero exit."""
    suite_data["samples"].append({"sample_id": "empty", "role": "mixed", "sex": "female"})
    suite = write_suite(suite_data)
    result = runner.invoke(app, ["build-truth", "--suite", str(suite)])
    assert result.exit_code == EXIT_SELECTION
    report = (_run_dir(suite) / "selection_report.md").read_text()
    assert "sample empty has no selected variants" in report
    assert not (_run_dir(suite) / "positive.truth.vcf.gz").exists()


def test_ref_mismatch_fails(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    """REF that does not match the FASTA -> non-zero exit with a clear message."""
    rule(suite_data, "synthetic_hom")["variants"].append(
        {
            "chrom": "chr2",
            "pos": 9800,
            "ref": "A",
            "alt": "G",
            "gene": "GENEC",
            "classification": "pathogenic",
            "rationale": "wrong REF on purpose",
        }
    )
    result = runner.invoke(app, ["build-truth", "--suite", str(write_suite(suite_data))])
    assert result.exit_code == EXIT_TRUTH
    assert "REF A does not match reference C at chr2:9800" in result.output


def test_manifest_from_other_config_is_refused(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path], tmp_path: Path
) -> None:
    out = tmp_path / "shared"
    assert runner.invoke(app, ["select", "--suite", str(SUITE), "--out", str(out)]).exit_code == 0
    suite_data["suite"]["seed"] = 1
    result = runner.invoke(
        app, ["select", "--suite", str(write_suite(suite_data)), "--out", str(out)]
    )
    assert result.exit_code == EXIT_INPUT
    assert "fresh output directory" in result.output


def test_validate_erepo_suite() -> None:
    result = runner.invoke(app, ["validate-suite", str(EREPO_SUITE)])
    assert result.exit_code == 0, result.output
    assert "fixture-erepo-suite 0.1 — 3 rules, 2 samples" in result.output


def test_erepo_reference_mismatch_fails(
    erepo_suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    """ERepo: a base stated in the HGVS differs from the FASTA -> exit 3, clear message."""
    erepo_suite_data["sources"]["clingen_erepo"]["file"] = str(FIXTURES / "erepo_badref.tsv")
    result = runner.invoke(app, ["select", "--suite", str(write_suite(erepo_suite_data))])
    assert result.exit_code == EXIT_TRUTH
    assert "eRepo HGVS reference bases differ from the reference FASTA" in result.output
    assert "NC_000002.12:g.9800A>G: stated reference A differs from the FASTA (C)" in result.output


def test_unsatisfied_erepo_rule_exit_code(
    erepo_suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    """ERepo: a required rule that cannot be satisfied -> exit 2 + report."""
    rule(erepo_suite_data, "erepo_plp")["per_gene"] = 5
    suite = write_suite(erepo_suite_data)
    result = runner.invoke(app, ["select", "--suite", str(suite)])
    assert result.exit_code == EXIT_SELECTION
    report = (_run_dir(suite) / "selection_report.md").read_text()
    assert "required rule erepo_plp is unsatisfied" in report
    assert "GENEX: requested 5, selected 2" in report


def test_simulation_commands_need_simulation_section(tmp_path: Path) -> None:
    for command in ("prepare-sim", "simulate"):
        result = runner.invoke(app, [command, "--suite", str(SUITE), "--out", str(tmp_path)])
        assert result.exit_code == EXIT_INPUT
        assert "no `simulation` section" in result.output


def test_run_logging_shows_progress_loggers_on_stderr(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with _run_logging(tmp_path, console=PROGRESS_LOGGERS):
        logging.getLogger("insilico_clinval.simulate.run").info("s1: lifting")
        logging.getLogger("insilico_clinval.variants.select").info("rule x: rejected y")
    err = capsys.readouterr().err
    assert "s1: lifting" in err and "rejected" not in err
    log_text = (tmp_path / LOG_NAME).read_text()
    assert "s1: lifting" in log_text and "rejected y" in log_text
    assert len(logging.getLogger("insilico_clinval").handlers) == 0
