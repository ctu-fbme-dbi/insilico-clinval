"""nf-core/sarek 3.10.0 on the synthetic fixture genome. Opt-in.

    CLINVAL_TEST_SAREK=1 CLINVAL_TEST_MASON_IMAGE=clinval-mason:2.4.0 \\
        pytest tests/integration/test_sarek.py

Needs Nextflow (host), Docker, network (sarek and its containers are pulled), bcftools,
samtools and Mason. Sarek builds the BWA index of the fixture FASTA itself and runs without base
recalibration (no known sites). Takes several minutes (sarek's containers are linux/amd64).
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pysam
import pytest
from typer.testing import CliRunner

from insilico_clinval.cli import app
from tests.conftest import EREPO_SUITE, FIXTURES, mason_runtime, requires_mason

pytestmark = [
    pytest.mark.integration,
    requires_mason,
    pytest.mark.skipif(
        os.environ.get("CLINVAL_TEST_SAREK") != "1" or shutil.which("nextflow") is None,
        reason="opt-in: CLINVAL_TEST_SAREK=1 and nextflow on PATH",
    ),
]


def test_sarek_calls_planted_snvs(tmp_path: Path) -> None:
    runner = CliRunner()
    out = tmp_path / "run"
    result = runner.invoke(app, ["build-truth", "--suite", str(EREPO_SUITE), "--out", str(out)])
    assert result.exit_code == 0, result.output
    args = ["simulate", "--suite", str(EREPO_SUITE), "--out", str(out)]
    runtime = mason_runtime(tmp_path)
    if runtime is not None:
        args += ["--runtime", str(runtime)]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    pipeline = FIXTURES / "pipeline_sarek_fixture.yaml"
    result = runner.invoke(
        app,
        ["run-sut", "--suite", str(EREPO_SUITE), "--out", str(out), "--pipeline", str(pipeline)],
    )
    assert result.exit_code == 0, result.output
    name = "sarek-3.10.0_fixture"
    for sample in ("positive", "negative"):
        secondary = json.loads((out / "sut" / name / f"{sample}.secondary.json").read_text())
        assert secondary["tool"] == "nf-core/sarek" and secondary["tool_version"] == "3.10.0"
        assert secondary["alignment"] is not None
        with pysam.VariantFile(str(out / secondary["raw_vcf"])) as vcf:
            called = {(r.chrom, r.pos, r.ref, r.alts[0]) for r in vcf if r.alts}
        expected = json.loads((out / f"{sample}.expected_outcomes.json").read_text())["outcomes"]
        snvs = []
        for outcome in expected:
            _, chrom, pos, ref, alt = outcome["variant_id"].split(":")
            if len(ref) == 1 and len(alt) == 1:
                snvs.append((chrom, int(pos), ref, alt))
        assert snvs
        missing = [v for v in snvs if v not in called]
        assert not missing, f"{sample}: planted SNVs not called by sarek: {missing}"
    manifest = json.loads((out / "run_manifest.json").read_text())
    tools = {t["name"]: t["version"] for t in manifest["tools"]}
    assert tools["nf-core/sarek"] == "3.10.0" and "nextflow" in tools
