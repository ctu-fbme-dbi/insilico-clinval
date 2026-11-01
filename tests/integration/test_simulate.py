"""Mason simulation of the eRepo fixture suite.

Needs bcftools, samtools and mason_simulator (on PATH, or a docker image named in
CLINVAL_TEST_MASON_IMAGE, see tests/conftest.py).
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pysam
import pytest
from typer.testing import CliRunner

from insilico_clinval.cli import app
from insilico_clinval.genome.reference import Reference
from insilico_clinval.simulate.check import variant_window
from tests.conftest import FIXTURES, mason_runtime, requires_mason

pytestmark = [pytest.mark.integration, requires_mason]

TRANS = ["curated:chr1:898:ATC:A", "curated:chr1:1499:T:TTGC"]
CIS = ["curated:chr2:600:GAG:TT", "curated:chr2:1200:T:A"]
HOM = "synthetic_novel:chr2:5100:G:T"
DETERMINISTIC = [
    "positive.R1.fastq.gz",
    "positive.R2.fastq.gz",
    "positive.true.bam",
    "positive.simulation.json",
    "positive.planting_check.json",
    "positive.sim/diploid.mason.vcf",
    "positive.sim/haploid.mason.vcf",
    "negative.true.bam",
]


def _simulate(suite: Path, out: Path, runtime: Path | None) -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["build-truth", "--suite", str(suite), "--out", str(out)])
    assert result.exit_code == 0, result.output
    args = ["simulate", "--suite", str(suite), "--out", str(out)]
    if runtime is not None:
        args += ["--runtime", str(runtime)]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    import yaml

    from tests.conftest import EREPO_SUITE, _absolutise

    base = tmp_path_factory.mktemp("sim")
    data = _absolutise(yaml.safe_load(EREPO_SUITE.read_text()))
    positive = next(s for s in data["samples"] if s["sample_id"] == "positive")
    positive["phase_sets"] = [
        {"variants": TRANS, "phase": "trans"},
        {"variants": CIS, "phase": "cis"},
    ]
    positive["overrides"] = [{"variant_id": HOM, "zygosity": "hom"}]
    data["output"] = {"directory": str(base / "default")}
    suite = base / "suite.yaml"
    suite.write_text(yaml.safe_dump(data, sort_keys=False))
    runtime = mason_runtime(base)
    a, b = base / "a", base / "b"
    _simulate(suite, a, runtime)
    _simulate(suite, b, runtime)
    return a, b


def test_mason_batches(tmp_path: Path) -> None:
    """simulation.mason_batch_bp: one Mason run per batch, merged as it ends, planting complete."""
    import yaml

    from tests.conftest import EREPO_SUITE, _absolutise

    data = _absolutise(yaml.safe_load(EREPO_SUITE.read_text()))
    data["simulation"]["mason_batch_bp"] = 1
    data["output"] = {"directory": str(tmp_path / "default")}
    suite = tmp_path / "suite.yaml"
    suite.write_text(yaml.safe_dump(data, sort_keys=False))
    out = tmp_path / "run"
    _simulate(suite, out, mason_runtime(tmp_path))
    spec = json.loads((out / "positive.simulation.json").read_text())
    names = [r["name"] for r in spec["runs"]]
    assert len(names) == len(set(names)) > 2
    assert all(len(r["segment_names"]) == 1 for r in spec["runs"])  # 1 bp: one segment each
    assert _check(out, "positive")
    sim_dir = out / "positive.sim"
    assert not list(sim_dir.glob("*.fq")) and not list(sim_dir.glob("*.true.bam"))
    assert sorted(p.name for p in sim_dir.glob("mason_*.log")) == sorted(
        f"mason_{n}.log" for n in names
    )
    with pysam.AlignmentFile(str(out / "positive.true.bam")) as bam:
        prefixes = {r.query_name.rsplit(".", 1)[0] for r in bam}
    assert prefixes == {f"positive.{n}" for n in names}


def _check(run: Path, sample: str) -> dict[str, dict[str, Any]]:
    data = json.loads((run / f"{sample}.planting_check.json").read_text())
    assert data["complete"] is True
    return {v["variant_id"]: v for v in data["variants"]}


def test_vaf_and_haplotypes(runs: tuple[Path, Path]) -> None:
    run, _ = runs
    for sample in ("positive", "negative"):
        for vid, v in _check(run, sample).items():
            if v["zygosity"] == "het":
                assert 0.3 <= v["vaf"] <= 0.7, (vid, v["vaf"])
            else:
                assert v["vaf"] >= 0.95, (vid, v["vaf"])
            for h in v["by_haplotype"]:
                if h["haplotype"] not in v["carrying_haplotypes"]:
                    assert h["alt"] <= 1, (vid, h)  # at most a coincidental sequencing error
                else:
                    assert h["alt"] >= 0.9 * h["reads"], (vid, h)


def test_phase_cis_trans_and_hom(runs: tuple[Path, Path]) -> None:
    check = _check(runs[0], "positive")
    trans = [check[v]["carrying_haplotypes"] for v in TRANS]
    cis = [check[v]["carrying_haplotypes"] for v in CIS]
    assert trans[0] != trans[1] and len(trans[0]) == len(trans[1]) == 1
    assert cis[0] == cis[1] and len(cis[0]) == 1
    assert check[HOM]["carrying_haplotypes"] == [1, 2]


def test_hemizygous_chrx_is_haploid_at_half_depth(runs: tuple[Path, Path]) -> None:
    run, _ = runs
    check = _check(run, "positive")
    hemi = [v["spanning_reads"] for v in check.values() if v["zygosity"] == "hemi"]
    het = [v["spanning_reads"] for v in check.values() if v["zygosity"] == "het"]
    assert hemi and het
    ratio = (sum(hemi) / len(hemi)) / (sum(het) / len(het))
    assert 0.3 <= ratio <= 0.75, ratio
    with pysam.AlignmentFile(str(run / "positive.true.bam")) as bam:
        haps = {r.get_tag("oH") for r in bam.fetch("chrX")}
        others = {r.get_tag("oH") for r in bam.fetch("chr1")}
    assert haps == {1}
    assert others == {1, 2}


def _net_length(read: pysam.AlignedSegment, start0: int, end0: int) -> int:
    """Inserted minus deleted bases of the alignment inside [start0, end0)."""
    net, ref_pos = 0, read.reference_start
    for op, n in read.cigartuples or []:
        if op == 1 and start0 <= ref_pos <= end0:
            net += n
        elif op == 2:
            net -= len(range(max(ref_pos, start0), min(ref_pos + n, end0)))
        if op in (0, 2, 3, 7, 8):
            ref_pos += n
    return net


def test_indel_lengths_are_exact(runs: tuple[Path, Path]) -> None:
    run, _ = runs
    check = _check(run, "positive")
    reference = Reference(FIXTURES / "ref.fa")
    with pysam.AlignmentFile(str(run / "positive.true.bam")) as bam:
        for vid, v in check.items():
            expected = len(v["alt"]) - len(v["ref"])
            start0, end0 = variant_window(reference, v["chrom"], v["pos"], v["ref"], v["alt"])
            carrying = [
                r
                for r in bam.fetch(v["chrom"], start0, end0)
                if r.get_tag("oH") in v["carrying_haplotypes"]
                and r.reference_start <= start0
                and (r.reference_end or 0) >= end0
            ]
            exact = sum(_net_length(r, start0, end0) == expected for r in carrying)
            assert carrying and exact >= 0.9 * len(carrying), (vid, exact, len(carrying))
    reference.close()


def test_same_seed_is_byte_identical(runs: tuple[Path, Path]) -> None:
    a, b = runs
    for name in DETERMINISTIC:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_manifest_records_simulation(runs: tuple[Path, Path]) -> None:
    manifest = json.loads((runs[0] / "run_manifest.json").read_text())
    tools = {t["name"]: t for t in manifest["tools"]}
    assert {"bcftools", "mason_simulator", "samtools"} <= set(tools)
    assert tools["mason_simulator"]["version"].startswith("mason_simulator version: 2.0.9")
    roles = {o["role"] for o in manifest["outputs"]}
    assert {"fastq_r1:positive", "true_bam:negative", "planting_check:positive"} <= roles
    assert "simulation_profile" in {i["role"] for i in manifest["inputs"]}
    spec = json.loads((runs[0] / "positive.simulation.json").read_text())
    assert [r["group"] for r in spec["runs"]] == ["diploid", "haploid"]
    assert spec["decomposed"] == ["curated:chr2:600:GAG:TT"]


def test_simulate_needs_truth(
    tmp_path: Path, write_suite: Callable[[dict[str, Any]], Path], erepo_suite_data: dict[str, Any]
) -> None:
    suite = write_suite(erepo_suite_data)
    result = CliRunner().invoke(
        app, ["simulate", "--suite", str(suite), "--out", str(tmp_path / "x")]
    )
    assert result.exit_code == 1
    assert "build-truth" in result.output


def test_recheck_planting_reproduces_the_check(runs: tuple[Path, Path], tmp_path: Path) -> None:
    """recheck-planting repeats the check from the true BAM and restores planting_complete."""
    import shutil

    run = tmp_path / "run"
    shutil.copytree(runs[0], run)
    original = (run / "positive.planting_check.json").read_bytes()
    spec_path = run / "positive.simulation.json"
    spec = json.loads(spec_path.read_text())
    spec["planting_complete"] = False
    spec_path.write_text(json.dumps(spec))
    suite = runs[0].parent / "suite.yaml"
    args = ["recheck-planting", "--suite", str(suite), "--out", str(run), "--sample", "positive"]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    assert (run / "positive.planting_check.json").read_bytes() == original
    assert json.loads(spec_path.read_text())["planting_complete"] is True
