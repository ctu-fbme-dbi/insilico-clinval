"""End-to-end build-truth with bcftools (genotypes, phasing, normalisation)."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pysam
import pytest
from typer.testing import CliRunner

from insilico_clinval.cli import EXIT_TRUTH, app
from insilico_clinval.provenance import sha256_file
from tests.conftest import SUITE, requires_bcftools, rule

pytestmark = [pytest.mark.integration, requires_bcftools]
runner = CliRunner()

DETERMINISTIC = [
    "variants.selected.json",
    "variants.selected.tsv",
    "selection_report.md",
    "positive.truth.vcf.gz",
    "positive.truth.vcf.gz.tbi",
    "positive.expected_outcomes.json",
    "positive.sample.json",
    "negative.truth.vcf.gz",
    "negative.truth.vcf.gz.tbi",
    "negative.expected_outcomes.json",
    "negative.sample.json",
]


def _build(suite: Path, out: Path) -> Path:
    result = runner.invoke(app, ["build-truth", "--suite", str(suite), "--out", str(out)])
    assert result.exit_code == 0, result.output
    return out


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _build(SUITE, tmp_path_factory.mktemp("bt") / "run")


def test_outputs_are_byte_identical(run_dir: Path, tmp_path: Path) -> None:
    """Same suite + seed -> byte-identical outputs (manifest/log excluded)."""
    again = _build(SUITE, tmp_path / "again")
    for name in DETERMINISTIC:
        assert sha256_file(run_dir / name) == sha256_file(again / name), name


def test_different_seed_changes_only_tie_breaks(
    run_dir: Path,
    suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    tmp_path: Path,
) -> None:
    suite_data["suite"]["seed"] = 7
    other = _build(write_suite(suite_data), tmp_path / "seed7")
    a = json.loads((run_dir / "variants.selected.json").read_text())
    b = json.loads((other / "variants.selected.json").read_text())
    assert [(r["rule_id"], r["selected"]) for r in a["rules"]] == [
        (r["rule_id"], r["selected"]) for r in b["rules"]
    ]
    assert len(a["variants"]) == len(b["variants"])


def test_truth_vcf_genotypes(run_dir: Path) -> None:
    """On the final (normalised, bgzipped) VCF."""
    with pysam.VariantFile(str(run_dir / "positive.truth.vcf.gz")) as vcf:
        recs = {(r.chrom, r.pos): r.samples["positive"] for r in vcf}
    assert recs[("chr2", 5100)]["GT"] == (1, 1)
    genea = [s for (c, p), s in recs.items() if c == "chr1" and p < 3000]
    assert len(genea) == 2
    assert {s["PS"] for s in genea} == {200}
    assert {s["GT"] for s in genea} == {(0, 1), (1, 0)}  # trans
    assert all(s.phased for s in recs.values())
    with pysam.VariantFile(str(run_dir / "negative.truth.vcf.gz")) as vcf:
        chrx = [r.samples["negative"]["GT"] for r in vcf if r.chrom == "chrX"]
    assert chrx == [(1,), (1,)]


def test_truth_vcf_header(run_dir: Path) -> None:
    """Each contig once; FORMAT only GT and PS; indexed."""
    with pysam.VariantFile(str(run_dir / "negative.truth.vcf.gz")) as vcf:
        assert list(vcf.header.contigs) == ["chr1", "chr2", "chrX", "chrY"]
        assert set(vcf.header.formats) == {"GT", "PS"}
        assert list(vcf.fetch("chrX"))  # tabix index works
    text = pysam.BGZFile(str(run_dir / "negative.truth.vcf.gz")).read().decode()
    assert text.count("##contig=<ID=chrX") == 1
    assert "bcftools_normCommand" not in text


def test_expected_outcomes_and_sample_spec(run_dir: Path) -> None:
    pos = json.loads((run_dir / "positive.expected_outcomes.json").read_text())
    neg = json.loads((run_dir / "negative.expected_outcomes.json").read_text())
    by_id = {o["variant_id"]: o for o in pos["outcomes"]}
    a1 = by_id["curated:chr1:200:T:G"]
    assert a1["expected_reportable"] is True
    assert a1["expected_classification"] == "pathogenic"
    assert a1["condition"] == {"id": "MONDO:0000101", "name": "Fixture disease A"}
    synth = by_id["synthetic_novel:chr2:5100:G:T"]
    assert synth["expected_acmg_criteria"] == ["PVS1", "PM2"]
    assert synth["evaluation_level"] == "interpretation"
    calling = [o for o in pos["outcomes"] if o["evaluation_level"] == "calling"]
    assert calling and all(o["expected_classification"] is None for o in calling)
    curated_neg = [o for o in neg["outcomes"] if o["variant_id"].startswith("curated:")]
    assert len(curated_neg) == 6
    assert all(not o["expected_reportable"] for o in curated_neg)
    assert all(o["reportable_basis"].startswith("negative control") for o in curated_neg)
    spec = json.loads((run_dir / "negative.sample.json").read_text())
    assert spec["truth_vcf"] == "negative.truth.vcf.gz"
    hemi = [v for v in spec["variants"] if v["zygosity"] == "hemi"]
    assert len(hemi) == 2


def test_manifest_provenance(run_dir: Path) -> None:
    """Manifest lists input SHA-256, outputs and the bcftools version."""
    m = json.loads((run_dir / "run_manifest.json").read_text())
    tools = {t["name"]: t["version"] for t in m["tools"]}
    assert tools["bcftools"].startswith("bcftools ")
    assert all(len(i["sha256"]) == 64 for i in m["inputs"])
    roles = {o["role"] for o in m["outputs"]}
    assert {"truth_vcf:positive", "expected_outcomes:negative", "sample_spec:positive"} <= roles
    names = [s["name"] for s in m["steps"]]
    assert names == [
        "select",
        "bcftools norm (positive)",
        "bcftools norm (negative)",
        "build-truth",
    ]
    assert all(s["exit_code"] == 0 for s in m["steps"])


def test_non_left_aligned_variant_is_an_error(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path], tmp_path: Path
) -> None:
    rule(suite_data, "synthetic_hom")["variants"].append(
        {
            "chrom": "chr2",
            "pos": 9503,
            "ref": "AA",
            "alt": "A",
            "gene": "GENEC",
            "classification": "pathogenic",
            "rationale": "deletion inside a homopolymer",
        }
    )
    out = tmp_path / "bad"
    result = runner.invoke(
        app, ["build-truth", "--suite", str(write_suite(suite_data)), "--out", str(out)]
    )
    assert result.exit_code == EXIT_TRUTH
    assert "bcftools norm changed the truth set" in result.output
    assert "9503" in result.output
