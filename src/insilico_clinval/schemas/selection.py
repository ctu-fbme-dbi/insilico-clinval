"""SelectionResult: persisted form of `variants.selected.json`."""

from pydantic import Field

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.variant import TestVariant


class Rejection(Contract):
    variant_id: str
    reason: str


class RuleSummary(Contract):
    rule_id: str
    required: bool
    samples: list[str]
    filter_counts: list[tuple[str, int]] = Field(
        description="(filter name, candidates remaining after it), in application order"
    )
    requested: dict[str, int] = Field(description="bucket (gene, bin, type, id, total) -> count")
    selected: dict[str, int]
    rejections: list[Rejection]
    missing: list[str]
    satisfied: bool


class SelectionResult(Contract):
    schema_version: int
    suite_name: str
    seed: int
    variants: list[TestVariant]
    rules: list[RuleSummary]
    warnings: list[str]
