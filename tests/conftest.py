"""Shared fixtures: the fixture suite and helpers to derive modified suites."""

from __future__ import annotations

import copy
import os
import shutil
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pytest
import yaml

from insilico_clinval.config import LoadedSuite, load_suite
from insilico_clinval.genome.reference import Reference

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SUITE = FIXTURES / "suite.yaml"
EREPO_SUITE = FIXTURES / "suite_erepo.yaml"

requires_bcftools = pytest.mark.skipif(shutil.which("bcftools") is None, reason="needs bcftools")

# Mason for integration tests: `mason_simulator` on PATH, or a docker image named in
# CLINVAL_TEST_MASON_IMAGE (e.g. clinval-mason:2.4.0 built from containers/mason).
MASON_IMAGE = os.environ.get("CLINVAL_TEST_MASON_IMAGE")
HAVE_MASON = shutil.which("mason_simulator") is not None or (
    MASON_IMAGE is not None and shutil.which("docker") is not None
)
requires_mason = pytest.mark.skipif(
    not HAVE_MASON or shutil.which("samtools") is None or shutil.which("bcftools") is None,
    reason="needs mason_simulator (PATH or CLINVAL_TEST_MASON_IMAGE), samtools and bcftools",
)


def mason_runtime(directory: Path) -> Path | None:
    """Runtime config for integration tests (None = everything local)."""
    if shutil.which("mason_simulator") is not None or MASON_IMAGE is None:
        return None
    path = directory / "runtime.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 1,
                "tools": {"mason_simulator": {"runner": "docker", "image": MASON_IMAGE}},
            }
        )
    )
    return path


def _absolutise(data: dict[str, Any]) -> dict[str, Any]:
    """Make every file path in a fixture-suite dict absolute (so it can live in tmp_path)."""
    data = copy.deepcopy(data)
    data["reference"]["fasta"] = str(FIXTURES / data["reference"]["fasta"])
    src = data["sources"]
    if "clingen_erepo" in src:
        erepo = src["clingen_erepo"]
        erepo["file"] = str(FIXTURES / erepo["file"])
        if erepo.get("held_out_reference"):
            ref = erepo["held_out_reference"]
            ref["file"] = str(FIXTURES / ref["file"])
    if data.get("background"):
        for key in ("vcf", "bed"):
            data["background"][key] = str(FIXTURES / data["background"][key])
    if "clinvar" in src:
        src["clinvar"]["vcf"] = str(FIXTURES / src["clinvar"]["vcf"])
    if "giab" in src:
        for key in ("truth_vcf", "confident_bed"):
            src["giab"][key] = str(FIXTURES / src["giab"][key])
    src["strata"] = {k: str(FIXTURES / v) for k, v in src.get("strata", {}).items()}
    if data.get("simulation"):
        data["simulation"]["profile"] = str(FIXTURES / data["simulation"]["profile"])
    return data


@pytest.fixture
def suite_data() -> dict[str, Any]:
    """The fixture suite as a mutable dict with absolute paths."""
    with SUITE.open() as fh:
        return _absolutise(yaml.safe_load(fh))


@pytest.fixture
def erepo_suite_data() -> dict[str, Any]:
    """The eRepo + background fixture suite as a mutable dict with absolute paths."""
    with EREPO_SUITE.open() as fh:
        return _absolutise(yaml.safe_load(fh))


@pytest.fixture
def write_suite(tmp_path: Path) -> Callable[[dict[str, Any]], Path]:
    counter = iter(range(1000))

    def _write(data: dict[str, Any]) -> Path:
        path = tmp_path / f"suite_{next(counter)}.yaml"
        data = copy.deepcopy(data)
        data["output"] = {"directory": str(tmp_path / f"run_{path.stem}")}
        path.write_text(yaml.safe_dump(data, sort_keys=False))
        return path

    return _write


@pytest.fixture
def loaded_suite() -> LoadedSuite:
    return load_suite(SUITE)


@pytest.fixture
def reference() -> Iterator[Reference]:
    ref = Reference(FIXTURES / "ref.fa")
    yield ref
    ref.close()


def rule(data: dict[str, Any], rule_id: str) -> dict[str, Any]:
    return next(r for r in data["selection"] if r["id"] == rule_id)
