"""Suite validation and JSON Schema export."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from insilico_clinval.cli import SCHEMA_MODELS, app
from insilico_clinval.config import SuiteError, load_suite
from tests.conftest import SUITE, rule

REPO = Path(__file__).resolve().parents[2]


def test_fixture_suite_is_valid() -> None:
    suite = load_suite(SUITE)
    assert suite.config.suite.seed == 20260928
    assert suite.config.selection[0].id == "clinvar_plp_panel"


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda d: d.update(unknown_key=1), "unknown_key"),
        (lambda d: rule(d, "clinvar_plp_panel").update(count=3), "exactly one of"),
        (lambda d: rule(d, "clinvar_plp_panel").update(samples=["nobody"]), "unknown samples"),
        (lambda d: rule(d, "difficult_cmrg").update(strata=["nope"]), "unknown stratum"),
        (
            lambda d: rule(d, "standard_giab").update(
                expected={"reportable": False, "level": "both"}
            ),
            "level must be calling",
        ),
        (lambda d: d["reference"].update(build="GRCh37"), "GRCh38"),
    ],
)
def test_invalid_suites(
    suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    mutate: Callable[[dict[str, Any]], None],
    message: str,
) -> None:
    mutate(suite_data)
    with pytest.raises(SuiteError, match=message):
        load_suite(write_suite(suite_data))


def test_duplicate_yaml_keys_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "dup.yaml"
    path.write_text("schema_version: 1\nschema_version: 1\n")
    with pytest.raises(SuiteError, match="duplicate key"):
        load_suite(path)


def test_exported_schemas_match_docs(tmp_path: Path) -> None:
    result = CliRunner().invoke(app, ["export-schemas", "--out", str(tmp_path)])
    assert result.exit_code == 0, result.output
    for name in SCHEMA_MODELS:
        exported = json.loads((tmp_path / f"{name}.schema.json").read_text())
        committed = json.loads((REPO / "docs" / "json_schema" / f"{name}.schema.json").read_text())
        assert exported == committed, f"docs/json_schema/{name}.schema.json is stale"


def test_example_suite_parses() -> None:
    suite = load_suite(REPO / "configs" / "suites" / "example.yaml")
    assert suite.resolve(suite.config.reference.fasta).name == "GRCh38_no_alt_analysis_set.fa"


def test_output_dir_collapses_parent_references(tmp_path: Path) -> None:
    suite = load_suite(SUITE)
    out = suite.output_dir(None)
    assert ".." not in out.parts
    assert out == SUITE.parent.parent.parent / "runs" / "fixture"
    assert suite.output_dir(tmp_path / "a" / ".." / "b") == tmp_path / "b"


def test_background_min_distance_defaults_to_constraints() -> None:
    config = load_suite(REPO / "configs" / "suites" / "example.yaml").config
    assert config.background is not None and config.background.min_distance_bp is None
    assert config.background_min_distance_bp() == config.constraints.min_distance_bp == 500
    own = config.model_copy(
        update={"background": config.background.model_copy(update={"min_distance_bp": 10})}
    )
    assert own.background_min_distance_bp() == 10
