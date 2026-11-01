"""scripts/background_diff.py on the fixture suite."""

from __future__ import annotations

import csv
import importlib.util
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any

from typer.testing import CliRunner

from insilico_clinval.cli import app

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "background_diff.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("background_diff", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["background_diff"] = module
    spec.loader.exec_module(module)
    return module


background_diff = _load()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_background_diff(
    tmp_path: Path,
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
) -> None:
    old = write_suite(erepo_suite_data)
    erepo_suite_data["background"]["restrict_to_bed"] = True
    new = write_suite(erepo_suite_data)
    run = tmp_path / "run"
    result = CliRunner().invoke(app, ["select", "--suite", str(old), "--out", str(run)])
    assert result.exit_code == 0, result.output

    outputs = []
    for name in ("a", "b"):
        out = tmp_path / name
        background_diff.main(
            ["--old", str(old), "--new", str(new), "--run", str(run), "--out", str(out)]
        )
        outputs.append({f: (out / f).read_bytes() for f in sorted(p.name for p in out.iterdir())})
    assert outputs[0] == outputs[1]  # deterministic
    out = tmp_path / "a"
    assert sorted(outputs[0]) == [
        "bed_sizes.tsv",
        "planted_context.tsv",
        "records.tsv",
        "skipped.tsv",
    ]

    sizes = {r["chrom"]: r for r in _rows(out / "bed_sizes.tsv")}
    assert all(r["old_bp"] == r["new_bp"] == r["both_bp"] for r in sizes.values())  # same BED
    records = _rows(out / "records.tsv")
    total = {
        b: sum(int(r["records"]) for r in records if (r["background"], r["stratum"]) == (b, "*"))
        for b in ("old", "new")
    }
    assert 0 < total["new"] <= total["old"]
    assert any(r["stratum"] == "cmrg" for r in records)
    skipped = {
        (r["background"], r["reason"]): int(r["records"]) for r in _rows(out / "skipped.tsv")
    }
    assert skipped.get(("new", "outside_bed"), 0) == total["old"] - total["new"]
    context = _rows(out / "planted_context.tsv")
    assert context and {"old_in_bed", "new_bg_within_150bp"} <= set(context[0])
