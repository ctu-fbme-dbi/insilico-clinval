"""TestVariant: one selected variant with full provenance."""

from pydantic import Field

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import Classification, VariantCategory, VariantType


class Source(Contract):
    database: str
    release: str
    record_id: str
    rsid: str | None = None


class Condition(Contract):
    name: str
    id: str | None = None


class ClinVarInfo(Contract):
    clnsig: list[str]
    review_status: str
    stars: int = Field(ge=0, le=4)
    conditions: list[Condition]


class ClinGenInfo(Contract):
    """ClinGen Evidence Repository (eRepo) assertion(s) for one allele.

    eRepo records of the same allele (e.g. curated for different conditions) are merged when their
    classifications fall in the same tier; ``erepo_id`` is the smallest Uuid, the others are listed
    in ``other_erepo_ids``.
    """

    vcep: str = Field(description="expert panel(s); several are joined with '; '")
    erepo_id: str
    other_erepo_ids: list[str] = Field(default_factory=list)
    caid: str | None = Field(default=None, description="ClinGen Allele Registry id")
    clinvar_variation_id: str | None = None
    erepo_release_first_seen: str = Field(
        description="held-out reference release if the allele is in it, else the current release"
    )
    classification: Classification = Field(description="most severe class of the merged records")
    met_criteria: list[str]
    not_met_criteria: list[str] = Field(default_factory=list)
    guideline: str | None = None
    conditions: list[Condition] = Field(default_factory=list)


class TestVariant(Contract):
    __test__ = False  # not a pytest test class

    variant_id: str
    rule_id: str
    chrom: str
    pos: int = Field(ge=1)
    ref: str
    alt: str
    variant_type: VariantType
    indel_length: int
    category: VariantCategory
    gene: str | None
    source: Source
    clinvar: ClinVarInfo | None
    clingen: ClinGenInfo | None = None
    held_out: bool = Field(
        default=False,
        description="absent from the held-out reference release (e.g. eRepo v2.2.0 for BIAS-2015)",
    )
    strata: list[str]
    selection_reason: str

    @staticmethod
    def make_id(category: VariantCategory, chrom: str, pos: int, ref: str, alt: str) -> str:
        return f"{category.value}:{chrom}:{pos}:{ref}:{alt}"
