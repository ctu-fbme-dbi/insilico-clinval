"""recheck-planting without a finished simulation. Needs bcftools (build-truth)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from insilico_clinval.cli import app
from tests.conftest import requires_bcftools

pytestmark = [pytest.mark.integration, requires_bcftools]


def test_recheck_planting_needs_simulation(
    tmp_path: Path, write_suite: Callable[[dict[str, Any]], Path], erepo_suite_data: dict[str, Any]
) -> None:
    suite = write_suite(erepo_suite_data)
    out = tmp_path / "x"
    runner = CliRunner()
    assert (
        runner.invoke(app, ["build-truth", "--suite", str(suite), "--out", str(out)]).exit_code == 0
    )
    result = runner.invoke(app, ["recheck-planting", "--suite", str(suite), "--out", str(out)])
    assert result.exit_code == 1
    assert "clinval simulate" in result.output
