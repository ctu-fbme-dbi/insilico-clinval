"""ExpectedOutcome (what the system under test should produce) and OutcomeRecord (what it did)."""

from typing import Literal

from pydantic import Field

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import Classification, EvaluationLevel, OutcomeClass
from insilico_clinval.schemas.variant import Condition


class ExpectedOutcome(Contract):
    variant_id: str
    expected_classification: Classification | None
    acceptable_classifications: list[Classification]
    expected_reportable: bool
    reportable_basis: str
    condition: Condition | None
    expected_acmg_criteria: list[str] = Field(default_factory=list)
    evaluation_level: EvaluationLevel
    primary_criterion: Literal["reportability", "classification"] = Field(
        default="reportability",
        description="2-tier reportability is the primary acceptance criterion; 5-tier is secondary",
    )
    notes: str = ""


class ExpectedOutcomeSet(Contract):
    """Persisted form of `<sample_id>.expected_outcomes.json`."""

    schema_version: int
    sample_id: str
    outcomes: list[ExpectedOutcome]


Evidence = dict[str, str | int | float | bool | list[str] | None]


class OutcomeRecord(Contract):
    """Comparison result for one planted variant, or one false report (``variant_id`` null)."""

    variant_id: str | None
    sample_id: str
    pipeline: str
    chrom: str
    pos: int
    ref: str
    alt: str
    evaluation_level: EvaluationLevel | None = Field(description="null for false reports")
    outcome_class: OutcomeClass
    stage: str | None = Field(description="last check passed, in pipeline order; null if none")
    expected_classification: Classification | None = None
    expected_reportable: bool | None = None
    observed_classification: Classification | None = None
    observed_reported: bool | None = None
    gene: str | None = Field(default=None, description="gene of the matched finding")
    evidence: Evidence = Field(default_factory=dict)


class OutcomeSet(Contract):
    """Persisted as `sut/<pipeline>/compare/<sample_id>.outcomes.json`."""

    schema_version: int
    sample_id: str
    pipeline: str
    parameters: dict[str, float | int | str]
    records: list[OutcomeRecord] = Field(description="one per planted variant, by variant_id")
    false_reports: list[OutcomeRecord] = Field(description="sorted by position")
