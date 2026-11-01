"""SuiteConfig: the YAML test-suite definition (see configs/suites/example.yaml).

File paths are plain strings here; they are resolved relative to the suite file by
`insilico_clinval.config`.
"""

from typing import Annotated, Literal

from pydantic import Field, model_validator

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import (
    Classification,
    EvaluationLevel,
    Phase,
    SampleRole,
    Sex,
    SimulationMode,
    VariantCategory,
    VariantType,
    Zygosity,
)
from insilico_clinval.schemas.variant import Condition


class SuiteInfo(Contract):
    name: str
    version: str
    seed: int


class Region(Contract):
    """1-based, closed interval."""

    chrom: str
    start: int = Field(ge=1)
    end: int = Field(ge=1)

    @model_validator(mode="after")
    def _ordered(self) -> "Region":
        if self.end < self.start:
            raise ValueError(f"region end < start: {self.chrom}:{self.start}-{self.end}")
        return self


class ReferenceConfig(Contract):
    build: Literal["GRCh38"]
    fasta: str
    sha256: str | None = None
    par_regions: list[Region] | None = Field(
        default=None,
        description="pseudo-autosomal regions; null = GRCh38 PAR1/PAR2 (see constants.py)",
    )


class ClinVarSource(Contract):
    vcf: str
    release: str


class GiabSource(Contract):
    truth_vcf: str
    confident_bed: str
    release: str
    name: str = "GIAB"


class HeldOutReference(Contract):
    """The held-out reference release; alleles absent from it are ``held_out``."""

    release: str = Field(description="e.g. v2.2.0 (the release BIAS-2015 was tuned on)")
    file: str = Field(
        description="eRepo export of that release (needs a Uuid or CAID column), or a GRCh38 VCF "
        "of its alleles (.vcf, .vcf.gz; matched by left-normalized allele)"
    )


class ClinGenErepoSource(Contract):
    file: str = Field(description="tab-delimited eRepo download")
    release: str
    held_out_reference: HeldOutReference | None = None


class Sources(Contract):
    clingen_erepo: ClinGenErepoSource | None = None
    clinvar: ClinVarSource | None = None
    giab: GiabSource | None = None
    strata: dict[str, str] = Field(default_factory=dict)


class BackgroundConfig(Contract):
    """Background genome (small-variant benchmark VCF + BED), e.g. GIAB HG002 v4.2.1 GRCh38."""

    name: str
    vcf: str = Field(description="bgzipped + tabix-indexed VCF")
    bed: str
    release: str
    sex: Sex = Field(description="sex of the background individual (HG002: male)")
    sample: str | None = Field(default=None, description="VCF sample column; null = the first")
    min_distance_bp: int | None = Field(
        default=None,
        ge=0,
        description="remove records closer than this to a planted variant; "
        "null = constraints.min_distance_bp",
    )
    # Record filters (defaults keep the v4.2.1 behaviour, so existing truth VCFs do not change).
    restrict_to_bed: bool = Field(
        default=False, description="skip records whose REF span is not inside the BED"
    )
    require_pass: bool = Field(
        default=False, description="skip records with a FILTER other than PASS or missing (.)"
    )
    max_length: int | None = Field(
        default=None, ge=1, description="skip records with a longer allele (event_length)"
    )
    spanning_deletion: Literal["skip", "ref"] = Field(
        default="skip",
        description="'*' alleles: skip the whole record, or read them as REF on their haplotype",
    )


class Constraints(Contract):
    min_distance_bp: int = Field(ge=0)
    allowed_variant_types: list[VariantType]
    max_indel_length: int | None = Field(default=None, ge=1)
    indel_ladder: list[tuple[int, int]] = Field(default_factory=list)

    @model_validator(mode="after")
    def _ladder(self) -> "Constraints":
        previous_end = 0
        for lo, hi in self.indel_ladder:
            if lo < 1 or hi < lo:
                raise ValueError(f"invalid indel ladder bin [{lo}, {hi}]")
            if lo <= previous_end:
                raise ValueError(f"indel ladder bins must be ascending and disjoint: [{lo}, {hi}]")
            previous_end = hi
        return self


class ExpectedSpec(Contract):
    classification: list[Classification] = Field(
        default_factory=list, description="acceptable classifications"
    )
    reportable: bool
    level: EvaluationLevel
    basis: str | None = None


class CompoundHet(Contract):
    """Pair selected variants of the same gene (in rank order) into phase sets."""

    phase: Phase


class _RuleBase(Contract):
    id: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    category: VariantCategory
    required: bool = False
    samples: list[str] = Field(min_length=1)
    expected: ExpectedSpec
    zygosity: Zygosity | None = None
    compound_het: CompoundHet | None = None


class ClinVarRule(_RuleBase):
    source: Literal["clinvar"]
    clnsig_allow: list[str] = Field(min_length=1)
    min_stars: int = Field(default=0, ge=0, le=4)
    clnvc: list[str] | None = None
    genes: list[str] | None = None
    conditions: list[str] | None = None
    region_bed: str | None = None
    strata: list[str] | None = None
    exclude_strata: list[str] | None = None
    variant_types: list[VariantType] | None = None
    # selection mode: exactly one of the following
    per_gene: int | None = Field(default=None, ge=1)
    count: int | None = Field(default=None, ge=1)
    ladder_per_bin: int | None = Field(default=None, ge=1)
    variation_ids: list[str] | None = None
    rsids: list[str] | None = None

    @model_validator(mode="after")
    def _mode(self) -> "ClinVarRule":
        modes = [
            m
            for m in ("per_gene", "count", "ladder_per_bin", "variation_ids", "rsids")
            if getattr(self, m) is not None
        ]
        if len(modes) != 1:
            raise ValueError(
                f"rule {self.id}: exactly one of per_gene, count, ladder_per_bin, "
                f"variation_ids, rsids is required (got {modes or 'none'})"
            )
        if self.per_gene is not None and not self.genes:
            raise ValueError(f"rule {self.id}: per_gene requires a gene list")
        if self.category not in (VariantCategory.CURATED, VariantCategory.DIFFICULT):
            raise ValueError(f"rule {self.id}: ClinVar rules must be curated or difficult")
        if (
            self.category == VariantCategory.DIFFICULT
            and not self.strata
            and self.ladder_per_bin is None
        ):
            raise ValueError(f"rule {self.id}: difficult rules need strata or ladder_per_bin")
        return self


class ClinGenRule(_RuleBase):
    source: Literal["clingen_erepo"]
    classification_allow: list[Classification] = Field(min_length=1)
    genes: list[str] | None = None
    vceps: list[str] | None = None
    conditions: list[str] | None = Field(default=None, description="MONDO ids")
    held_out: Literal["include", "exclude", "only"] = "include"
    region_bed: str | None = None
    strata: list[str] | None = None
    exclude_strata: list[str] | None = None
    variant_types: list[VariantType] | None = None
    # selection mode: exactly one of the following
    per_gene: int | None = Field(default=None, ge=1)
    count: int | None = Field(default=None, ge=1)
    ladder_per_bin: int | None = Field(default=None, ge=1)
    erepo_ids: list[str] | None = None

    @model_validator(mode="after")
    def _mode(self) -> "ClinGenRule":
        modes = [
            m
            for m in ("per_gene", "count", "ladder_per_bin", "erepo_ids")
            if getattr(self, m) is not None
        ]
        if len(modes) != 1:
            raise ValueError(
                f"rule {self.id}: exactly one of per_gene, count, ladder_per_bin, erepo_ids is "
                f"required (got {modes or 'none'})"
            )
        if self.per_gene is not None and not self.genes:
            raise ValueError(f"rule {self.id}: per_gene requires a gene list")
        if self.category not in (VariantCategory.CURATED, VariantCategory.DIFFICULT):
            raise ValueError(f"rule {self.id}: eRepo rules must be curated or difficult")
        if (
            self.category == VariantCategory.DIFFICULT
            and not self.strata
            and self.ladder_per_bin is None
        ):
            raise ValueError(f"rule {self.id}: difficult rules need strata or ladder_per_bin")
        return self


class GiabRule(_RuleBase):
    source: Literal["giab"]
    region_bed: str | None = None
    strata: list[str] | None = None
    exclude_strata: list[str] | None = None
    variant_types: list[VariantType] | None = None
    counts: dict[VariantType, int] | None = None
    count: int | None = Field(default=None, ge=1)
    ladder_per_bin: int | None = Field(default=None, ge=1)

    @model_validator(mode="after")
    def _mode(self) -> "GiabRule":
        modes = [m for m in ("counts", "count", "ladder_per_bin") if getattr(self, m) is not None]
        if len(modes) != 1:
            raise ValueError(
                f"rule {self.id}: exactly one of counts, count, ladder_per_bin is required"
            )
        if self.counts is not None and any(n < 1 for n in self.counts.values()):
            raise ValueError(f"rule {self.id}: counts must be positive")
        if self.category not in (VariantCategory.STANDARD, VariantCategory.DIFFICULT):
            raise ValueError(f"rule {self.id}: GIAB rules must be standard or difficult")
        if (
            self.category == VariantCategory.DIFFICULT
            and not self.strata
            and self.ladder_per_bin is None
        ):
            raise ValueError(f"rule {self.id}: difficult rules need strata or ladder_per_bin")
        if self.expected.level != EvaluationLevel.CALLING:
            raise ValueError(
                f"rule {self.id}: GIAB variants carry no classification; level must be calling"
            )
        if self.compound_het is not None:
            raise ValueError(f"rule {self.id}: compound_het needs gene annotation (ClinVar)")
        return self


class ExplicitVariant(Contract):
    chrom: str
    pos: int = Field(ge=1)
    ref: str = Field(pattern=r"^[ACGTacgt]+$")
    alt: str = Field(pattern=r"^[ACGTacgt]+$")
    gene: str | None = None
    classification: Classification | None = None
    condition: Condition | None = None
    acmg_criteria: list[str] = Field(default_factory=list)
    rationale: str


class ExplicitRule(_RuleBase):
    source: Literal["explicit"]
    variants: list[ExplicitVariant]

    @model_validator(mode="after")
    def _explicit(self) -> "ExplicitRule":
        if self.category != VariantCategory.SYNTHETIC_NOVEL:
            raise ValueError(f"rule {self.id}: explicit rules must be synthetic_novel")
        if self.expected.level != EvaluationLevel.CALLING:
            for v in self.variants:
                if v.classification is None:
                    raise ValueError(
                        f"rule {self.id}: variant {v.chrom}:{v.pos} needs an intended "
                        f"classification for level {self.expected.level.value}"
                    )
        return self


Rule = Annotated[ClinGenRule | ClinVarRule | GiabRule | ExplicitRule, Field(discriminator="source")]
# Rules whose variants carry a gene annotation and region/strata filters.
AnnotatedRule = ClinGenRule | ClinVarRule


class VariantOverride(Contract):
    variant_id: str
    zygosity: Zygosity


class PhaseSetSpec(Contract):
    variants: list[str] = Field(min_length=2)
    phase: Phase

    @model_validator(mode="after")
    def _trans_pair(self) -> "PhaseSetSpec":
        if self.phase == Phase.TRANS and len(self.variants) != 2:
            raise ValueError("a trans phase set must contain exactly two variants")
        return self


class SampleConfig(Contract):
    sample_id: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    role: SampleRole
    sex: Sex
    default_zygosity: Literal[Zygosity.HET, Zygosity.HOM] = Zygosity.HET
    overrides: list[VariantOverride] = Field(default_factory=list)
    phase_sets: list[PhaseSetSpec] = Field(default_factory=list)
    background: bool | None = Field(
        default=None, description="merge the suite background; null = yes if one is configured"
    )


class OutputConfig(Contract):
    directory: str


class SimulationConfig(Contract):
    """Read simulation: profile file and which part of the genome is simulated."""

    profile: str = Field(description="simulation profile YAML (configs/profiles/)")
    mode: SimulationMode = SimulationMode.REGION
    flank_bp: int = Field(
        default=1500, ge=0, description="region mode: bases simulated on each side of a variant"
    )
    extra_regions_bed: str | None = Field(
        default=None, description="region mode: additional regions, e.g. paralogous loci"
    )
    genome_contigs: list[str] | None = Field(
        default=None, description="genome mode: contigs to simulate; null = chr1-22, chrX, chrY"
    )
    min_alt_reads: int = Field(
        default=1, ge=1, description="planting check: alt reads required on each carrying haplotype"
    )
    mason_batch_bp: int | None = Field(
        default=None,
        ge=1,
        description=(
            "split each ploidy group into Mason runs of consecutive whole segments of at most this "
            "many bp (a longer segment runs alone); caps Mason's temporary disk use; changes the "
            "reads; null = one Mason run per group"
        ),
    )

    @model_validator(mode="after")
    def _mode(self) -> "SimulationConfig":
        if self.mode == SimulationMode.GENOME and self.extra_regions_bed is not None:
            raise ValueError("simulation.extra_regions_bed applies to region mode only")
        if self.mode == SimulationMode.REGION and self.genome_contigs is not None:
            raise ValueError("simulation.genome_contigs applies to genome mode only")
        return self


class SuiteConfig(Contract):
    schema_version: Literal[1]
    suite: SuiteInfo
    reference: ReferenceConfig
    sources: Sources
    constraints: Constraints
    selection: list[Rule] = Field(min_length=1)
    background: BackgroundConfig | None = None
    samples: list[SampleConfig] = Field(min_length=1)
    simulation: SimulationConfig | None = None
    output: OutputConfig

    def uses_background(self, sample: SampleConfig) -> bool:
        if self.background is None:
            return False
        return sample.background is not False

    def background_min_distance_bp(self) -> int:
        """Gap kept between background records and planted variants."""
        if self.background is not None and self.background.min_distance_bp is not None:
            return self.background.min_distance_bp
        return self.constraints.min_distance_bp

    @model_validator(mode="after")
    def _cross_references(self) -> "SuiteConfig":
        errors: list[str] = []
        rule_ids = [r.id for r in self.selection]
        sample_ids = [s.sample_id for s in self.samples]
        for kind, ids in (("rule", rule_ids), ("sample", sample_ids)):
            duplicates = sorted({i for i in ids if ids.count(i) > 1})
            if duplicates:
                errors.append(f"duplicate {kind} ids: {', '.join(duplicates)}")
        for sample in self.samples:
            if sample.background and self.background is None:
                errors.append(f"sample {sample.sample_id}: background requested but not configured")
            if self.uses_background(sample):
                assert self.background is not None
                if sample.sex != self.background.sex:
                    errors.append(
                        f"sample {sample.sample_id}: sex {sample.sex.value} differs from the "
                        f"{self.background.name} background ({self.background.sex.value})"
                    )
        erepo = self.sources.clingen_erepo
        for rule in self.selection:
            unknown = [s for s in rule.samples if s not in sample_ids]
            if unknown:
                errors.append(f"rule {rule.id}: unknown samples {unknown}")
            if isinstance(rule, ClinGenRule):
                if erepo is None:
                    errors.append(f"rule {rule.id}: sources.clingen_erepo is not configured")
                elif rule.held_out != "include" and erepo.held_out_reference is None:
                    errors.append(
                        f"rule {rule.id}: held_out: {rule.held_out} needs "
                        "sources.clingen_erepo.held_out_reference"
                    )
            if isinstance(rule, ClinVarRule) and self.sources.clinvar is None:
                errors.append(f"rule {rule.id}: sources.clinvar is not configured")
            if isinstance(rule, GiabRule) and self.sources.giab is None:
                errors.append(f"rule {rule.id}: sources.giab is not configured")
            if isinstance(rule, (ClinGenRule, ClinVarRule, GiabRule)):
                for name in (rule.strata or []) + (rule.exclude_strata or []):
                    if name not in self.sources.strata:
                        errors.append(f"rule {rule.id}: unknown stratum {name!r}")
                if rule.ladder_per_bin is not None and not self.constraints.indel_ladder:
                    errors.append(f"rule {rule.id}: ladder_per_bin needs constraints.indel_ladder")
            if rule.expected.level != EvaluationLevel.CALLING and not rule.expected.classification:
                errors.append(
                    f"rule {rule.id}: expected.classification is required for level "
                    f"{rule.expected.level.value}"
                )
        if errors:
            raise ValueError("; ".join(errors))
        return self
