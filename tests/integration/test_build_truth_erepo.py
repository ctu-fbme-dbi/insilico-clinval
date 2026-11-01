"""build-truth with eRepo variants and the HG002-like background."""

from __future__ import annotations

import csv
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pysam
import pytest
from typer.testing import CliRunner

from insilico_clinval._runner import Toolbox
from insilico_clinval.cli import app
from insilico_clinval.genome.alleles import ref_span0
from insilico_clinval.provenance import sha256_file
from insilico_clinval.truth.normalize import normalize_vcf
from tests.conftest import EREPO_SUITE, FIXTURES, requires_bcftools

pytestmark = [pytest.mark.integration, requires_bcftools]
runner = CliRunner()

SAMPLES = ("positive", "negative")
DETERMINISTIC = [
    "variants.selected.json",
    "variants.selected.tsv",
    "selection_report.md",
    *(
        f"{s}.{suffix}"
        for s in SAMPLES
        for suffix in (
            "truth.vcf.gz",
            "truth.vcf.gz.tbi",
            "expected_outcomes.json",
            "sample.json",
            "background_removed.tsv",
        )
    ),
]


def _build(suite: Path, out: Path) -> Path:
    result = runner.invoke(app, ["build-truth", "--suite", str(suite), "--out", str(out)])
    assert result.exit_code == 0, result.output
    return out


def _records(path: Path, sample: str) -> list[tuple[str, int, str, tuple[str, ...], Any, bool]]:
    with pysam.VariantFile(str(path)) as vcf:
        return [
            (
                r.chrom,
                r.pos,
                str(r.ref),
                tuple(r.alts or ()),
                r.samples[sample]["GT"],
                "TV_BG" in r.info,
            )
            for r in vcf
        ]


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _build(EREPO_SUITE, tmp_path_factory.mktemp("erepo") / "run")


def test_outputs_are_byte_identical(run_dir: Path, tmp_path: Path) -> None:
    """Same suite + seed -> byte-identical outputs, including the background merge."""
    again = _build(EREPO_SUITE, tmp_path / "again")
    for name in DETERMINISTIC:
        assert sha256_file(run_dir / name) == sha256_file(again / name), name


def test_different_seed_changes_only_tie_breaks(
    run_dir: Path,
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    tmp_path: Path,
) -> None:
    erepo_suite_data["suite"]["seed"] = 11
    other = _build(write_suite(erepo_suite_data), tmp_path / "seed11")
    a = json.loads((run_dir / "variants.selected.json").read_text())
    b = json.loads((other / "variants.selected.json").read_text())
    assert [(r["rule_id"], r["selected"]) for r in a["rules"]] == [
        (r["rule_id"], r["selected"]) for r in b["rules"]
    ]
    # background records not near planted variants are identical apart from haplotype choice
    bg_a = {
        (c, p, r, al)
        for c, p, r, al, _, bg in _records(run_dir / "negative.truth.vcf.gz", "negative")
        if bg
    }
    bg_b = {
        (c, p, r, al)
        for c, p, r, al, _, bg in _records(other / "negative.truth.vcf.gz", "negative")
        if bg
    }
    assert len(bg_a ^ bg_b) <= 4


def test_background_merge(run_dir: Path) -> None:
    """Background records near planted variants are removed and logged."""
    pos_recs = _records(run_dir / "positive.truth.vcf.gz", "positive")
    neg_recs = _records(run_dir / "negative.truth.vcf.gz", "negative")
    # chr2:5100 G>A (background) overlaps the synthetic chr2:5100 G>T of the positive sample
    at_5100 = [r for r in pos_recs if r[:2] == ("chr2", 5100)]
    assert [(r[3], r[5]) for r in at_5100] == [(("T",), False)]
    assert any(r[5] and r[:2] == ("chr2", 5100) for r in neg_recs)
    with (run_dir / "positive.background_removed.tsv").open() as fh:
        removed = list(csv.DictReader(fh, delimiter="\t"))
    assert {
        "chrom": "chr2",
        "pos": "5100",
        "ref": "G",
        "alt": "A",
        "gt": "0|1",
        "near_variant_id": "synthetic_novel:chr2:5100:G:T",
        "gap_bp": "-1",
    } in removed
    for sample, recs in (("positive", pos_recs), ("negative", neg_recs)):
        planted = [(c, *ref_span0(p, r)) for c, p, r, _, _, bg in recs if not bg]
        for c, p, r, _, _, bg in recs:
            if not bg:
                continue
            start, end = ref_span0(p, r)
            gaps = [max(start - e, s - end) for pc, s, e in planted if pc == c]
            assert all(g >= 300 for g in gaps), (sample, c, p)  # min_distance_bp
    spec = json.loads((run_dir / "positive.sample.json").read_text())
    bg = spec["background"]
    assert bg["name"] == "GIAB_HG002" and bg["release"] == "fixture-v4.2.1"
    assert bg["removed"] == len(removed) and bg["removed_log"] == "positive.background_removed.tsv"
    # chrX (no background truth in v4.2.1 either) and chr2:9499+ lie outside the BED
    assert all(
        v.split(":")[1] == "chrX" or int(v.split(":")[2]) > 9000 for v in bg["planted_outside_bed"]
    )
    assert any(v.split(":")[1] == "chrX" for v in bg["planted_outside_bed"])


def test_background_genotypes(run_dir: Path) -> None:
    """Unphased background het -> phased; multi-allelic 1|2 split; chrX hemizygous."""
    recs = _records(run_dir / "negative.truth.vcf.gz", "negative")
    bg = {(c, p, al): gt for c, p, _, al, gt, is_bg in recs if is_bg}
    assert bg[("chr1", 7000, ("G",))] == (1, 1)
    assert set(bg) >= {("chr1", 7600, ("G",)), ("chr1", 7600, ("A",))}
    assert {bg[("chr1", 7600, ("G",))], bg[("chr1", 7600, ("A",))]} == {(1, 0), (0, 1)}
    assert not any(p in (8000, 8500) and c == "chr1" and al == ("T",) for c, p, al in bg)  # 0/0
    assert ("chr2", 8500, ("C",)) not in bg  # ./.
    with pysam.VariantFile(str(run_dir / "negative.truth.vcf.gz")) as vcf:
        for r in vcf:
            call = r.samples["negative"]
            if r.chrom == "chrX":
                assert call["GT"] == (1,)
            else:
                assert call.phased
    # chr1:6500 is unphased 0/1 in the background; it is phased in every sample that keeps it
    pos = {
        (c, p): gt
        for c, p, _, _, gt, is_bg in _records(run_dir / "positive.truth.vcf.gz", "positive")
        if is_bg
    }
    assert pos[("chr1", 6500)] in {(0, 1), (1, 0)}


def test_header_and_normalisation(run_dir: Path, tmp_path: Path) -> None:
    """Contigs once, FORMAT GT/PS only, sorted and unchanged by bcftools norm."""
    path = run_dir / "positive.truth.vcf.gz"
    with pysam.VariantFile(str(path)) as vcf:
        assert list(vcf.header.contigs) == ["chr1", "chr2", "chrX", "chrY"]
        assert set(vcf.header.formats) == {"GT", "PS"}
        assert {"TV_BG", "TV_HELD_OUT", "TV_VCEP"} <= set(vcf.header.info)
    text = pysam.BGZFile(str(path)).read().decode()
    assert text.count("##contig=<ID=chr1") == 1
    assert "##clinval_background=<Name=GIAB_HG002,Release=fixture-v4.2.1,Sample=HG002>" in text
    plain = tmp_path / "positive.vcf"
    plain.write_text(text)
    again = normalize_vcf(
        plain,
        tmp_path / "again.vcf",
        FIXTURES / "ref.fa",
        step="t",
        manifest=None,
        toolbox=Toolbox(),
    )
    with pysam.VariantFile(str(plain)) as a, pysam.VariantFile(str(again)) as b:
        before = [str(r) for r in a]
        after = [str(r) for r in b]
    assert before == after  # all records, background included, are already normalised
    keys = [
        (("chr1", "chr2", "chrX", "chrY").index(line.split("\t")[0]), int(line.split("\t")[1]))
        for line in before
    ]
    assert keys == sorted(keys)


def test_expected_outcomes_for_erepo(run_dir: Path) -> None:
    sel = {
        v["variant_id"]: v
        for v in json.loads((run_dir / "variants.selected.json").read_text())["variants"]
    }
    pos = json.loads((run_dir / "positive.expected_outcomes.json").read_text())["outcomes"]
    neg = json.loads((run_dir / "negative.expected_outcomes.json").read_text())["outcomes"]
    for o in pos + neg:
        v = sel[o["variant_id"]]
        if v["clingen"] is None:
            continue
        assert o["expected_classification"] == v["clingen"]["classification"]
        assert o["expected_acmg_criteria"] == v["clingen"]["met_criteria"]
        assert o["condition"]["id"].startswith("MONDO:")
        assert v["clingen"]["erepo_id"] in o["notes"]
        assert ("held out: absent from eRepo v2.2.0" in o["notes"]) is v["held_out"]
    assert all(o["expected_reportable"] for o in pos)
    assert not any(o["expected_reportable"] for o in neg)
    assert all(o["reportable_basis"].startswith("negative control") for o in neg)


def test_manifest_provenance(run_dir: Path) -> None:
    """Input SHA-256 (eRepo, held-out reference, background) and bcftools version."""
    m = json.loads((run_dir / "run_manifest.json").read_text())
    inputs = {i["role"]: i["sha256"] for i in m["inputs"]}
    assert inputs["clingen_erepo"] == sha256_file(FIXTURES / "erepo_mini.tsv")
    assert inputs["clingen_erepo_held_out"] == sha256_file(FIXTURES / "erepo_v2.2.0_mini.tsv")
    assert inputs["background_vcf"] == sha256_file(FIXTURES / "hg002_mini.vcf.gz")
    assert {"background_vcf_index", "background_bed"} <= set(inputs)
    dbs = {d["name"]: d["release"] for d in m["databases"]}
    assert dbs == {
        "ClinGen eRepo": "fixture-2026-09",
        "ClinGen eRepo (held-out reference)": "v2.2.0",
        "GIAB_HG002": "fixture-v4.2.1",
    }
    assert {t["name"] for t in m["tools"]} == {"bcftools"}
    assert {"background_removed:positive", "background_removed:negative"} <= {
        o["role"] for o in m["outputs"]
    }
    assert any("background records removed" in w for w in m["warnings"])
