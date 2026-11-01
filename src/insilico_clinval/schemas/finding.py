"""Normalised outputs of the system under test: SecondaryResult and CanonicalFinding."""

from pydantic import Field

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import Classification


class SecondaryResult(Contract):
    """Persisted as `sut/<pipeline>/<sample_id>.secondary.json`; paths relative to the run dir."""

    schema_version: int
    sample_id: str
    pipeline: str
    adapter: str
    tool: str
    tool_version: str
    alignment: str | None
    raw_vcf: str
    filtered_vcf: str | None


class CanonicalFinding(Contract):
    sample_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    classification: Classification | None
    reported: bool
    acmg_criteria: list[str] = Field(default_factory=list)
    tool: str
    tool_version: str


class FindingSet(Contract):
    """Persisted as `sut/<pipeline>/<sample_id>.findings.json` (findings sorted by position)."""

    schema_version: int
    sample_id: str
    pipeline: str
    adapter: str
    tool: str
    tool_version: str
    input_vcf: str = Field(description="tertiary input VCF, relative to the run directory")
    findings: list[CanonicalFinding]
