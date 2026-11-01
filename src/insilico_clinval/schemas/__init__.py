"""
Pydantic data contracts shared across modules. Source of truth: docs/schemas.md.
"""

from insilico_clinval.schemas.enums import (
    Classification,
    EvaluationLevel,
    Haplotype,
    OutcomeClass,
    Phase,
    SampleRole,
    Sex,
    VariantCategory,
    VariantType,
    Zygosity,
)
from insilico_clinval.schemas.manifest import (
    ConfigRef,
    DatabaseRef,
    FileRef,
    FrameworkInfo,
    RunManifest,
    StepRecord,
    ToolRef,
)
from insilico_clinval.schemas.outcome import (
    ExpectedOutcome,
    ExpectedOutcomeSet,
    OutcomeRecord,
    OutcomeSet,
)
from insilico_clinval.schemas.sample import BackgroundInfo, SampleSpec, SampleVariant
from insilico_clinval.schemas.selection import RuleSummary, SelectionResult
from insilico_clinval.schemas.suite import (
    ClinGenRule,
    ClinVarRule,
    ExplicitRule,
    ExplicitVariant,
    GiabRule,
    Rule,
    SuiteConfig,
)
from insilico_clinval.schemas.variant import (
    ClinGenInfo,
    ClinVarInfo,
    Condition,
    Source,
    TestVariant,
)

SCHEMA_VERSION = 1

__all__ = [
    "SCHEMA_VERSION",
    "BackgroundInfo",
    "Classification",
    "ClinGenInfo",
    "ClinGenRule",
    "ClinVarInfo",
    "ClinVarRule",
    "Condition",
    "ConfigRef",
    "DatabaseRef",
    "EvaluationLevel",
    "ExpectedOutcome",
    "ExpectedOutcomeSet",
    "ExplicitRule",
    "ExplicitVariant",
    "FileRef",
    "FrameworkInfo",
    "GiabRule",
    "Haplotype",
    "OutcomeClass",
    "OutcomeRecord",
    "OutcomeSet",
    "Phase",
    "Rule",
    "RuleSummary",
    "RunManifest",
    "SampleRole",
    "SampleSpec",
    "SampleVariant",
    "SelectionResult",
    "Sex",
    "Source",
    "StepRecord",
    "SuiteConfig",
    "TestVariant",
    "ToolRef",
    "VariantCategory",
    "VariantType",
    "Zygosity",
]
