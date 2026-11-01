"""CLNSIG allow-list, review stars, conditions, genes."""

from __future__ import annotations

import pytest

from insilico_clinval.schemas.enums import Classification
from insilico_clinval.variants.clinvar import (
    clnsig_tokens,
    expected_class,
    parse_conditions,
    parse_genes,
    passes_allow_list,
)
from insilico_clinval.variants.stars import REVIEW_STARS, UnknownReviewStatus, review_stars

PLP = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}


@pytest.mark.parametrize(
    ("status", "stars"),
    [
        ("practice_guideline", 4),
        ("reviewed_by_expert_panel", 3),
        ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
        ("criteria_provided,_single_submitter", 1),
        ("criteria_provided,_conflicting_classifications", 1),
        ("no_assertion_criteria_provided", 0),
        ("no_classification_provided", 0),
    ],
)
def test_review_stars(status: str, stars: int) -> None:
    assert review_stars(status) == stars


def test_every_mapped_status_is_in_range() -> None:
    assert set(REVIEW_STARS.values()) <= {0, 1, 2, 3, 4}


def test_unknown_review_status_is_an_error() -> None:
    with pytest.raises(UnknownReviewStatus):
        review_stars("made_up_status")


def test_clnsig_tokens() -> None:
    assert clnsig_tokens("Pathogenic|Uncertain_significance") == [
        "pathogenic",
        "uncertain_significance",
    ]
    assert clnsig_tokens("Pathogenic,_low_penetrance") == ["pathogenic", "low_penetrance"]
    assert clnsig_tokens("Pathogenic/Likely_pathogenic") == ["pathogenic/likely_pathogenic"]


def test_strict_allow_list_rejects_mixed() -> None:
    assert passes_allow_list(["pathogenic"], PLP)
    assert passes_allow_list(["pathogenic/likely_pathogenic"], PLP)
    assert not passes_allow_list(clnsig_tokens("Pathogenic|Uncertain_significance"), PLP)
    assert not passes_allow_list(clnsig_tokens("Pathogenic|risk_factor"), PLP)
    assert not passes_allow_list(["conflicting_classifications_of_pathogenicity"], PLP)
    assert not passes_allow_list([], PLP)


def test_expected_class() -> None:
    assert expected_class(["pathogenic/likely_pathogenic"]) == Classification.PATHOGENIC
    assert expected_class(["likely_benign"]) == Classification.LIKELY_BENIGN
    assert expected_class(["benign/likely_benign"]) == Classification.BENIGN
    assert expected_class(["risk_factor"]) is None


def test_parse_conditions_prefers_mondo_and_skips_placeholders() -> None:
    conds = parse_conditions(
        "MONDO:MONDO:0000101,OMIM:600101|MedGen:CN169374|OMIM:600200,MedGen:C1",
        "Disease_A|not_specified|Disease_B",
    )
    assert [(c.name, c.id) for c, _ in conds] == [
        ("Disease A", "MONDO:0000101"),
        ("Disease B", "OMIM:600200"),
    ]
    assert conds[0][1] == {"MONDO:0000101", "OMIM:600101"}


def test_parse_genes() -> None:
    assert parse_genes("GENEB:9002|GENEZ:9099") == ["GENEB", "GENEZ"]
    assert parse_genes(None) == []
