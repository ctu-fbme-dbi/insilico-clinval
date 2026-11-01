"""Variant selection on the fixture suite (counts, rules, tie-breaks, bins)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from insilico_clinval.config import LoadedSuite, load_suite
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.selection import SelectionResult
from insilico_clinval.variants.report import selected_tsv, selection_report
from insilico_clinval.variants.select import (
    Deduplicator,
    select_variants,
    selection_failures,
    variants_per_sample,
)
from tests.conftest import rule


@pytest.fixture
def result(loaded_suite: LoadedSuite, reference: Reference) -> SelectionResult:
    return select_variants(loaded_suite, reference)


def _summary(result: SelectionResult, rule_id: str):  # type: ignore[no-untyped-def]
    return next(r for r in result.rules if r.rule_id == rule_id)


def _record_ids(result: SelectionResult, rule_id: str) -> set[str]:
    return {v.source.record_id for v in result.variants if v.rule_id == rule_id}


def test_fixture_selection_counts(result: SelectionResult) -> None:
    assert len(result.variants) == 25
    per_sample = variants_per_sample(result)
    assert len(per_sample["positive"]) == 19
    assert len(per_sample["negative"]) == 11
    assert selection_failures(result, ["positive", "negative"]) == []


def test_practice_guideline_wins_ranking(result: SelectionResult) -> None:
    """The 4-star record is selected and is the top-ranked GENEA candidate."""
    a1 = next(v for v in result.variants if v.source.record_id == "100001")
    assert a1.clinvar is not None and a1.clinvar.stars == 4
    assert "rank 1/" in a1.selection_reason
    assert a1.source.rsid == "rs1001"
    # 1-star GENEA pathogenic record never beats higher-star candidates
    assert "100006" not in _record_ids(result, "clinvar_plp_panel")


def test_mixed_clnsig_rejected(result: SelectionResult) -> None:
    """Pathogenic|Uncertain_significance is not selected by a P/LP allow-list."""
    assert all(v.source.record_id != "100005" for v in result.variants)
    counts = dict(_summary(result, "clinvar_plp_panel").filter_counts)
    assert counts["clnsig_allow"] < counts["max_indel_length"]


def test_too_close_candidate_replaced(result: SelectionResult) -> None:
    """100002 (3 stars, 250 bp from 100001) is rejected; the next candidate is taken."""
    summary = _summary(result, "clinvar_plp_panel")
    rejected = {r.variant_id: r.reason for r in summary.rejections}
    assert "min_distance_bp=500" in rejected["curated:chr1:450:A:C"]
    genea = [v for v in result.variants if v.rule_id == "clinvar_plp_panel" and v.gene == "GENEA"]
    assert len(genea) == 2
    # exactly one of the overlapping pair 100003/100004 is in the suite
    ids = {v.source.record_id for v in result.variants}
    assert len(ids & {"100003", "100004"}) == 1


def test_overlap_rejected_across_rules(result: SelectionResult) -> None:
    ladder = _summary(result, "indel_ladder")
    reasons = [r.reason for r in ladder.rejections]
    assert any("already selected" in r or "overlaps" in r for r in reasons)


def test_deduplicator() -> None:
    d = Deduplicator(min_distance_bp=10)
    assert d.conflict("chr1", 100, "ACG", "A") is None
    from insilico_clinval.schemas.enums import VariantCategory, VariantType
    from insilico_clinval.schemas.variant import Source, TestVariant

    tv = TestVariant(
        variant_id="curated:chr1:100:ACG:A",
        rule_id="r",
        chrom="chr1",
        pos=100,
        ref="ACG",
        alt="A",
        variant_type=VariantType.DELETION,
        indel_length=-2,
        category=VariantCategory.CURATED,
        gene=None,
        source=Source(database="x", release="1", record_id="1"),
        clinvar=None,
        strata=[],
        selection_reason="",
    )
    d.accept(tv)
    assert "same variant" in (d.conflict("chr1", 100, "ACG", "A") or "")
    assert "overlaps" in (d.conflict("chr1", 101, "C", "T") or "")
    assert "min_distance_bp=10" in (d.conflict("chr1", 108, "C", "T") or "")  # gap 5
    assert d.conflict("chr1", 113, "C", "T") is None  # gap 10
    assert d.conflict("chr2", 100, "A", "T") is None


def test_indel_ladder_shortfall(result: SelectionResult) -> None:
    """Requested count per bin, shortfall reported for 21-50 bp."""
    ladder = _summary(result, "indel_ladder")
    assert ladder.selected == {"1-5bp": 2, "6-10bp": 2, "11-20bp": 2, "21-50bp": 1}
    assert not ladder.satisfied
    assert any(m.startswith("21-50bp: requested 2, selected 1") for m in ladder.missing)
    assert "100028" not in _record_ids(result, "indel_ladder")  # 60 bp > max_indel_length


def test_provenance_fields(result: SelectionResult) -> None:
    """Every TestVariant carries source database, release and record id."""
    for v in result.variants:
        assert v.source.database and v.source.release and v.source.record_id
    clinvar = [v for v in result.variants if v.source.database == "ClinVar"]
    assert {v.source.release for v in clinvar} == {"2026-09-01"}


def test_strata_and_contig_mapping(result: SelectionResult) -> None:
    overlapping = next(v for v in result.variants if v.source.record_id in {"100003", "100004"})
    assert overlapping.chrom == "chr1"  # ClinVar contig "1" mapped to the reference name
    assert overlapping.strata == ["cmrg"]
    assert any("contig_not_in_reference" in w for w in result.warnings)


def test_male_and_giab_filters(result: SelectionResult) -> None:
    giab = [v for v in result.variants if v.source.database == "GIAB"]
    assert {v.pos for v in giab if v.rule_id == "standard_giab"} <= {6400, 7000, 7100, 7800}
    assert all(v.strata == ["cmrg"] for v in giab if v.rule_id == "difficult_cmrg")
    assert all(v.pos not in (8500, 9000, 9200) for v in giab)


def test_different_seed_changes_only_tie_breaks(
    suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
    result: SelectionResult,
) -> None:
    """Counts and non-tied choices are seed-independent."""
    suite_data["suite"]["seed"] = 1
    other = select_variants(load_suite(write_suite(suite_data)), reference)
    assert [(r.rule_id, r.selected) for r in other.rules] == [
        (r.rule_id, r.selected) for r in result.rules
    ]
    for record_id in ("100001", "100011", "100012", "100031", "100032"):
        assert any(v.source.record_id == record_id for v in other.variants)


def test_same_seed_same_outputs(loaded_suite: LoadedSuite, reference: Reference) -> None:
    a = select_variants(loaded_suite, reference)
    b = select_variants(loaded_suite, reference)
    assert a == b
    assert selected_tsv(a) == selected_tsv(b)
    assert selection_report(loaded_suite, a) == selection_report(loaded_suite, b)


def test_explicit_ids_and_rsids(
    suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path], reference: Reference
) -> None:
    r = rule(suite_data, "clinvar_plp_panel")
    for key in ("genes", "per_gene", "compound_het"):
        r.pop(key)
    r["variation_ids"] = ["100011", "999999"]
    res = select_variants(load_suite(write_suite(suite_data)), reference)
    summary = _summary(res, "clinvar_plp_panel")
    assert summary.selected == {"VariationID 100011": 1, "VariationID 999999": 0}
    assert not summary.satisfied

    r.pop("variation_ids")
    r["rsids"] = ["rs1001"]
    res = select_variants(load_suite(write_suite(suite_data)), reference)
    assert _record_ids(res, "clinvar_plp_panel") == {"100001"}
