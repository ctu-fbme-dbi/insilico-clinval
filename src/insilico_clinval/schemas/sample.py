"""SampleSpec: one simulated genome."""

from pydantic import Field

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import Haplotype, SampleRole, Sex, Zygosity


class SampleVariant(Contract):
    variant_id: str
    zygosity: Zygosity
    haplotype: Haplotype
    phase_set: int | None = Field(
        default=None, description="VCF PS value (position of the first member); null if none"
    )


class BackgroundInfo(Contract):
    """Background genome merged into the truth VCF (e.g. GIAB HG002 v4.2.1)."""

    name: str
    release: str
    vcf: str = Field(description="as configured in the suite")
    bed: str = Field(description="as configured in the suite")
    records: int = Field(description="background records written to the truth VCF (before split)")
    removed: int = Field(description="records removed within min_distance_bp of planted variants")
    min_distance_bp: int | None = Field(
        default=None, description="gap used for the removal (null in runs built before it was set)"
    )
    removed_log: str = Field(description="TSV of removed records, relative to the run directory")
    planted_outside_bed: list[str] = Field(
        description="planted variant ids outside the background high-confidence BED"
    )
    skipped: dict[str, int] = Field(
        default_factory=dict, description="background records skipped, by reason"
    )
    spanning_deletion_to_ref: int = Field(
        default=0, description="'*' alleles read as REF (spanning_deletion: ref)"
    )


class SampleSpec(Contract):
    schema_version: int
    sample_id: str
    role: SampleRole
    sex: Sex
    background: BackgroundInfo | None = None
    variants: list[SampleVariant]
    truth_vcf: str = Field(description="path relative to the run directory")
    expected_outcomes: str = Field(description="path relative to the run directory")
