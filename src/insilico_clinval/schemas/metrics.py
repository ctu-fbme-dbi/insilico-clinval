"""Metrics, acceptance criteria and regression diff, see docs/schemas.md.

The acceptance configuration is separate from the suite, so defining or changing the criteria
does not invalidate a run directory.
"""

from typing import Literal

from pydantic import Field, model_validator

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.fidelity import Proportion

Subset = dict[str, str]

PROPORTION_METRICS = (
    "report_sensitivity", "report_specificity", "class_exact", "class_acceptable", "class_tier",
    "criterion_recall", "criterion_precision", "planted_calling", "calling_recall",
    "calling_precision",
)  # fmt: skip
COUNT_METRICS = ("stage_counts", "false_reports", "false_positive_calls")


# -- acceptance configuration (YAML input) ------------------------------------------------------


class Criterion(Contract):
    metric: str
    subset: Subset = Field(default_factory=dict, description="e.g. {held_out: 'true'}")
    min_point: float | None = Field(default=None, ge=0, le=1)
    min_lower: float | None = Field(default=None, ge=0, le=1, description="lower Wilson bound")
    max_count: int | None = Field(default=None, ge=0)

    @model_validator(mode="after")
    def _check(self) -> "Criterion":
        if self.metric not in (*PROPORTION_METRICS, *COUNT_METRICS):
            raise ValueError(f"unknown metric {self.metric!r}")
        bounds = (self.min_point, self.min_lower, self.max_count)
        if all(b is None for b in bounds):
            raise ValueError(f"criterion {self.metric}: set min_point, min_lower or max_count")
        is_count = self.metric in COUNT_METRICS
        if is_count and (self.min_point is not None or self.min_lower is not None):
            raise ValueError(f"criterion {self.metric}: a count takes max_count only")
        if not is_count and self.max_count is not None:
            raise ValueError(f"criterion {self.metric}: a proportion takes min_point/min_lower")
        return self


class AcceptanceConfig(Contract):
    schema_version: Literal[1] = 1
    name: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    criteria: list[Criterion] = Field(min_length=1)


# -- results ------------------------------------------------------------------------------------


class MetricRow(Contract):
    metric: str
    subset: Subset
    value: Proportion


class CountRow(Contract):
    metric: str
    subset: Subset
    count: int


class CriterionResult(Contract):
    criterion: Criterion
    status: Literal["pass", "fail", "not_evaluable"]
    observed: str = Field(description="the value the criterion was judged on")


class MetricsReport(Contract):
    """Persisted as `metrics.json`."""

    schema_version: int
    pipeline: str
    runs: list[str] = Field(description="run directories, as given")
    proportions: list[MetricRow]
    counts: list[CountRow]
    acceptance: str | None = Field(default=None, description="acceptance config name")
    criteria: list[CriterionResult] = Field(default_factory=list)


ChangeKind = Literal[
    "new_failure", "recovered", "class_change", "stage_change", "false_report_added",
    "false_report_removed",
]  # fmt: skip


class Change(Contract):
    kind: ChangeKind
    run: str
    sample_id: str
    variant_id: str | None
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str | None
    baseline: str = Field(description="baseline outcome (class / observed class / reported)")
    candidate: str
    noise: bool = Field(description="also seen between the baseline and the noise run")


class CallingDelta(Contract):
    run: str
    sample_id: str
    contig_group: str
    counts: dict[str, int] = Field(description="truth TP/FN, query FP of baseline and candidate")


class Concordance(Contract):
    run: str
    sample_id: str
    contig_group: str
    identical: int
    only_baseline: int
    only_candidate: int


class DiffReport(Contract):
    """Persisted as `diff.json`."""

    schema_version: int
    baseline: str
    candidate: str
    noise: str | None
    runs: list[str]
    changes: list[Change]
    calling: list[CallingDelta]
    concordance: list[Concordance] = Field(
        description="PASS calls of both raw VCFs; empty when both use the same secondary outputs"
    )
    verdict: Literal["no_regression", "regression"]
    acceptance: str | None = None
    criteria: list[CriterionResult] = Field(default_factory=list)
