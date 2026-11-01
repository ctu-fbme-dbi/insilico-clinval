"""ClinVar review status (CLNREVSTAT) to review stars.

Mapping follows the ClinVar documentation "Review status in ClinVar"
(https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/). Values are the VCF spellings
(spaces replaced by underscores). ``criteria_provided,_conflicting_classifications`` has one
star in ClinVar (see docs/schemas.md). Older releases spell it
``..._conflicting_interpretations``. An unknown value is an error, never silently 0.
"""

REVIEW_STARS: dict[str, int] = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_classification_for_the_single_variant": 0,
    "no_interpretation_for_the_single_variant": 0,
    "no_classifications_from_unflagged_records": 0,
}


class UnknownReviewStatus(ValueError):
    pass


def review_stars(review_status: str) -> int:
    key = review_status.strip().lower()
    try:
        return REVIEW_STARS[key]
    except KeyError:
        raise UnknownReviewStatus(
            f"unknown ClinVar review status {review_status!r}; extend "
            "insilico_clinval.variants.stars.REVIEW_STARS"
        ) from None
