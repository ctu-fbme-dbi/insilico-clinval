"""Failure-stage attribution for one planted variant (pure logic, no I/O).

The checks run in pipeline order and the first failing check sets the outcome class. Checks
beyond the variant's evaluation level are not scored: a ``calling`` variant that is called with
the right genotype is ``correct``. A pipeline without tertiary analysis is scored as if every
variant had the ``calling`` level.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from insilico_clinval.schemas.enums import (
    Classification,
    EvaluationLevel,
    OutcomeClass,
    PlantingStatus,
)
from insilico_clinval.schemas.outcome import Evidence, ExpectedOutcome

# Check names (OutcomeRecord.stage = last check passed), in pipeline order.
CHECKS = (
    "simulated", "covered", "mapped", "called", "passed", "genotyped",
    "annotated", "classified", "reported",
)  # fmt: skip
CALLING_CHECKS = 6


@dataclass(frozen=True)
class VariantEvidence:
    planting_status: PlantingStatus
    true_alt_reads: int
    min_alt_reads: int
    aligned_alt_reads: int | None  # None: the pipeline exposes no alignment
    happy_decision: str | None  # TRUTH BD: TP, FN, UNK; None: allele not in hap.py output
    happy_kind: str | None  # TRUTH BK: gm, am, lm, .
    nonpass_calls: list[str] = field(default_factory=list)  # FILTER of non-PASS records
    tertiary: bool = True  # the pipeline has a tertiary stage
    finding_found: bool = False
    observed_classification: Classification | None = None
    observed_reported: bool | None = None


@dataclass(frozen=True)
class Attribution:
    outcome: OutcomeClass
    stage: str | None


def _failures(
    expected: ExpectedOutcome, ev: VariantEvidence, mismapped_fraction: float
) -> list[OutcomeClass | None]:
    """Failure class of each check in CHECKS order (None = passed)."""
    called = ev.happy_decision == "TP"
    # the read count only explains a missed allele: aligners clip reads around long indels, so a
    # called allele can have few aligned reads that carry it
    allele_called = called or ev.happy_kind == "am"
    fails: list[OutcomeClass | None] = [
        OutcomeClass.NOT_SIMULATED if ev.planting_status != PlantingStatus.OK else None,
        OutcomeClass.NOT_COVERED if ev.true_alt_reads < ev.min_alt_reads else None,
        OutcomeClass.MISMAPPED
        if not allele_called
        and ev.aligned_alt_reads is not None
        and ev.aligned_alt_reads < mismapped_fraction * ev.true_alt_reads
        else None,
    ]
    if called:
        fails += [None, None, None]
    elif ev.happy_kind == "am":  # a PASS call of the allele with the wrong genotype
        fails += [None, None, OutcomeClass.GENOTYPE_MISMATCH]
    elif ev.nonpass_calls:
        fails += [None, OutcomeClass.FILTERED, None]
    else:
        fails += [OutcomeClass.NOT_CALLED, None, None]
    acceptable = set(expected.acceptable_classifications)
    fails += [
        OutcomeClass.NOT_ANNOTATED
        if not ev.finding_found or ev.observed_classification is None
        else None,
        OutcomeClass.MISCLASSIFIED if ev.observed_classification not in acceptable else None,
        OutcomeClass.REPORT_MISMATCH
        if ev.observed_reported != expected.expected_reportable
        else None,
    ]
    return fails


def attribute(
    expected: ExpectedOutcome, ev: VariantEvidence, *, mismapped_fraction: float
) -> Attribution:
    scored = len(CHECKS)
    if expected.evaluation_level == EvaluationLevel.CALLING or not ev.tertiary:
        scored = CALLING_CHECKS
    fails = _failures(expected, ev, mismapped_fraction)[:scored]
    for i, failure in enumerate(fails):
        if failure is not None:
            return Attribution(failure, CHECKS[i - 1] if i else None)
    return Attribution(OutcomeClass.CORRECT, CHECKS[scored - 1])


def evidence_dict(ev: VariantEvidence, **extra: str | int | float | bool | None) -> Evidence:
    out: Evidence = {
        "planting_status": ev.planting_status.value,
        "true_alt_reads": ev.true_alt_reads,
        "aligned_alt_reads": ev.aligned_alt_reads,
        "happy_decision": ev.happy_decision,
        "happy_kind": ev.happy_kind,
        "nonpass_filters": list(ev.nonpass_calls),
        "finding_found": ev.finding_found,
    }
    out.update(extra)
    return out
