"""The README quickstart (examples/quickstart/run_quickstart.sh) end to end. Opt-in, needs the
images of configs/runtime/docker-local-build.yaml:

CLINVAL_TEST_MASON_IMAGE=clinval-mason:2.4.0 CLINVAL_TEST_HAPPY_IMAGE=clinval-happy:0.3.15 \
    pytest tests/integration/test_quickstart.py
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from tests.conftest import requires_mason

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "examples" / "quickstart" / "run_quickstart.sh"

pytestmark = [
    pytest.mark.integration,
    requires_mason,
    pytest.mark.skipif(
        os.environ.get("CLINVAL_TEST_HAPPY_IMAGE") is None or shutil.which("docker") is None,
        reason="needs hap.py image",
    ),
]


def test_quickstart(tmp_path: Path) -> None:
    out = tmp_path / "quickstart"
    env = {**os.environ, "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}"}
    result = subprocess.run(
        ["bash", str(SCRIPT), str(out)], env=env, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr[-3000:]
    with open(out / "sut" / "toy" / "compare" / "outcomes.tsv", newline="") as fh:
        rows = {r["variant_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    # the two failures the README describes
    delins, novel = rows["curated:chr2:600:GAG:TT"], rows["synthetic_novel:chr2:5100:G:T"]
    assert (delins["outcome_class"], delins["stage"]) == ("not_called", "mapped")
    assert (novel["outcome_class"], novel["stage"]) == ("misclassified", "annotated")
    assert sum(r["outcome_class"] == "correct" for r in rows.values()) == len(rows) - 2
    diff = json.loads((out / "diff" / "diff.json").read_text())
    assert diff["verdict"] == "regression"
    assert "Regression diff" in (out / "report.md").read_text()
