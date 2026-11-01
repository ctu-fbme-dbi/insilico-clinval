"""Expected-outcome metadata derived from the selection rule of each variant."""

from __future__ import annotations

from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import Classification, EvaluationLevel, SampleRole
from insilico_clinval.schemas.outcome import ExpectedOutcome, ExpectedOutcomeSet
from insilico_clinval.schemas.suite import (
    ClinGenRule,
    ClinVarRule,
    ExplicitRule,
    ExplicitVariant,
    SampleConfig,
    SuiteConfig,
)
from insilico_clinval.schemas.variant import Condition, TestVariant
from insilico_clinval.truth import TruthError
from insilico_clinval.truth.genotype import Genotype
from insilico_clinval.variants.clinvar import expected_class

_BENIGN = {Classification.BENIGN, Classification.LIKELY_BENIGN}


def _explicit_entry(rule: ExplicitRule, v: TestVariant) -> ExplicitVariant:
    for entry in rule.variants:
        if (entry.chrom, entry.pos, entry.ref.upper(), entry.alt.upper()) == (
            v.chrom,
            v.pos,
            v.ref,
            v.alt,
        ):
            return entry
    raise TruthError(f"{v.variant_id}: not found in explicit rule {rule.id}")


def _default_basis(
    rule_id: str, level: EvaluationLevel, reportable: bool, cls: Classification | None
) -> str:
    if level == EvaluationLevel.CALLING:
        return f"calling-only variant (rule {rule_id}); interpretation not scored"
    if reportable:
        return f"{cls.value if cls else 'classified'} variant in scope of rule {rule_id}"
    if cls is not None and cls in _BENIGN:
        return f"negative control: {cls.value} variant must not be reported (rule {rule_id})"
    return f"must not be reported (rule {rule_id})"


def _clinvar_condition(rule: ClinVarRule, v: TestVariant) -> Condition | None:
    if v.clinvar is None or not v.clinvar.conditions:
        return None
    conditions = v.clinvar.conditions
    if rule.conditions:
        wanted = set(rule.conditions)
        for c in conditions:
            if c.id in wanted:
                return c
    return next((c for c in conditions if c.id is not None), conditions[0])


def _erepo_condition(rule: ClinGenRule, v: TestVariant) -> Condition | None:
    assert v.clingen is not None
    conditions = v.clingen.conditions
    if rule.conditions:
        wanted = set(rule.conditions)
        for c in conditions:
            if c.id in wanted:
                return c
    return conditions[0] if conditions else None


def _erepo_note(cfg: SuiteConfig, v: TestVariant) -> str:
    assert v.clingen is not None
    src = cfg.sources.clingen_erepo
    ref = src.held_out_reference if src is not None else None
    note = f"ClinGen eRepo {v.clingen.erepo_id} ({v.clingen.vcep})"
    if ref is None:
        return f"{note}; held-out status not determined (no held_out_reference)"
    if v.held_out:
        return f"{note}; held out: absent from eRepo {ref.release}"
    return f"{note}; present in eRepo {ref.release}"


def expected_outcome(cfg: SuiteConfig, g: Genotype) -> ExpectedOutcome:
    v = g.variant
    rule = next(r for r in cfg.selection if r.id == v.rule_id)
    spec = rule.expected
    acceptable = list(spec.classification)
    cls: Classification | None = None
    condition: Condition | None = None
    criteria: list[str] = []
    notes = list(g.notes)
    if isinstance(rule, ClinGenRule):
        assert v.clingen is not None
        condition = _erepo_condition(rule, v)
        criteria = list(v.clingen.met_criteria)
        if spec.level != EvaluationLevel.CALLING:
            cls = v.clingen.classification
        notes.append(_erepo_note(cfg, v))
    elif isinstance(rule, ClinVarRule):
        condition = _clinvar_condition(rule, v)
        if spec.level != EvaluationLevel.CALLING:
            assert v.clinvar is not None
            cls = expected_class(v.clinvar.clnsig)
    elif isinstance(rule, ExplicitRule):
        entry = _explicit_entry(rule, v)
        condition = entry.condition
        criteria = list(entry.acmg_criteria)
        if spec.level != EvaluationLevel.CALLING:
            cls = entry.classification
        notes.append(f"rationale: {entry.rationale}")
    if spec.level != EvaluationLevel.CALLING:
        if cls is None:
            raise TruthError(f"{v.variant_id}: no expected classification can be derived")
        if cls not in acceptable:
            raise TruthError(
                f"{v.variant_id}: expected classification {cls.value} is not in the acceptable "
                f"classifications of rule {rule.id} {[c.value for c in acceptable]}"
            )
    else:
        acceptable = []
    return ExpectedOutcome(
        variant_id=v.variant_id,
        expected_classification=cls,
        acceptable_classifications=acceptable,
        expected_reportable=spec.reportable,
        reportable_basis=spec.basis or _default_basis(rule.id, spec.level, spec.reportable, cls),
        condition=condition,
        expected_acmg_criteria=criteria,
        evaluation_level=spec.level,
        notes="; ".join(notes),
    )


def expected_outcomes(
    cfg: SuiteConfig, sample: SampleConfig, genotypes: list[Genotype]
) -> tuple[ExpectedOutcomeSet, list[str]]:
    outcomes = sorted((expected_outcome(cfg, g) for g in genotypes), key=lambda o: o.variant_id)
    warnings: list[str] = []
    n_reportable = sum(o.expected_reportable for o in outcomes)
    if sample.role == SampleRole.NEGATIVE and n_reportable:
        warnings.append(
            f"negative sample {sample.sample_id} has {n_reportable} expected reportable findings"
        )
    if sample.role == SampleRole.POSITIVE and not n_reportable:
        warnings.append(f"positive sample {sample.sample_id} has no expected reportable findings")
    result = ExpectedOutcomeSet(
        schema_version=SCHEMA_VERSION, sample_id=sample.sample_id, outcomes=outcomes
    )
    return result, warnings
