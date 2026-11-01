"""Fidelity experiment: configuration and result contracts, see docs/schemas.md.

The configuration is separate from the suite: it compares secondary-analysis outputs of a
simulated and a real genome. Paths are strings, resolved relative to the configuration file.
"""

from typing import Literal

from pydantic import Field, model_validator

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import PlantingStatus, Sex, Zygosity

# -- configuration (YAML input) -------------------------------------------------------------


class FidelityDataset(Contract):
    id: str = Field(pattern=r"^[A-Za-z0-9_-]+$", description="file names and table key")
    label: str = Field(description="shown in tables and figures")
    alignment: str = Field(description="BAM or CRAM of the secondary pipeline")
    vcf: str = Field(description="VCF of the secondary pipeline (bgzipped, indexed)")
    truth_vcf: str = Field(description="truth of this genome (bgzipped, indexed)")


class HappySettings(Contract):
    engine: Literal["vcfeval", "xcmp"] = "vcfeval"
    threads: int = Field(default=4, ge=1)
    pass_only: bool = Field(default=True, description="evaluate PASS calls only")
    vcfeval_template: str | None = Field(
        default=None, description="RTG SDF of the reference; null = hap.py builds one"
    )
    gender: Sex | None = Field(
        default=None, description="hap.py --gender (HG002: male); null = hap.py guesses"
    )


class FidelityConfig(Contract):
    schema_version: Literal[1] = 1
    name: str = Field(pattern=r"^[A-Za-z0-9_.-]+$")
    reference: str = Field(description="FASTA (with .fai) used by the pipeline")
    confident_bed: str = Field(description="regions evaluated by hap.py (GIAB benchmark BED)")
    qc_regions: str | None = Field(
        default=None, description="regions of samtools/bcftools stats; null = confident_bed"
    )
    stratifications: dict[str, str] = Field(
        default_factory=dict, description="stratum name -> BED (plain or .gz)"
    )
    datasets: list[FidelityDataset] = Field(min_length=2)
    baseline: str = Field(description="dataset id the others are compared with (real data)")
    recall_gap: float | None = Field(
        default=None, gt=0, lt=1, description="flag strata whose |recall gap| exceeds this"
    )
    happy: HappySettings = Field(default_factory=HappySettings)
    planting: list[str] = Field(
        default_factory=list, description="run directories with *.planting_check.json"
    )
    threads: int = Field(default=4, ge=1, description="samtools / bcftools threads")

    @model_validator(mode="after")
    def _ids(self) -> "FidelityConfig":
        ids = [d.id for d in self.datasets]
        if len(set(ids)) != len(ids):
            raise ValueError(f"dataset ids must be unique: {ids}")
        if self.baseline not in ids:
            raise ValueError(f"baseline {self.baseline!r} is not a dataset id ({', '.join(ids)})")
        for name in self.stratifications:
            if not name or "\t" in name or name == "*":
                raise ValueError(f"invalid stratum name {name!r}")
        return self


# -- results (`fidelity.json`) --------------------------------------------------------------


class Proportion(Contract):
    """k of n with a 95 % Wilson score interval; value and bounds are null when n = 0."""

    k: int
    n: int
    value: float | None
    lower: float | None
    upper: float | None


class ReadQc(Contract):
    """From `samtools stats` in the QC regions (error and indel rates include true variants)."""

    dataset: str
    reads_total: int
    reads_mapped: float = Field(description="fraction of reads")
    reads_properly_paired: float
    reads_mq0: float
    insert_size_mean: float
    insert_size_sd: float
    error_rate: float = Field(description="mismatches / bases mapped (cigar)")
    indel_rate: float = Field(description="indel events / bases mapped (cigar)")
    coverage_mean: float
    coverage_ge20: float = Field(description="fraction of region bases with depth >= 20")
    fold80: float | None = Field(description="mean / 20th percentile of depth; null if it is 0")


class CallsetQc(Contract):
    """From `bcftools stats` in the QC regions."""

    dataset: str
    snvs: int
    indels: int
    ti_tv: float | None
    het_hom: float | None


class CallingMetric(Contract):
    """One hap.py row: variant type x stratum (subtype `*`, genotype `*`)."""

    dataset: str
    type: Literal["SNP", "INDEL"]
    stratum: str = Field(description="`*` = all confident regions")
    subtype: str = Field(description="hap.py subtype, `*` = all")
    recall: Proportion
    precision: Proportion
    f1: float | None


class RecallGap(Contract):
    dataset: str
    type: Literal["SNP", "INDEL"]
    stratum: str
    recall: float | None
    baseline_recall: float | None
    gap: float | None = Field(description="recall - baseline recall")
    flagged: bool


class PlantingRow(Contract):
    run: str = Field(description="run directory name")
    sample_id: str
    variant_id: str
    zygosity: Zygosity
    expected_vaf: float
    spanning_reads: int
    alt_reads: int
    vaf: float | None
    status: PlantingStatus


class ZygositySummary(Contract):
    zygosity: Zygosity
    variants: int
    vaf_median: float | None
    vaf_min: float | None
    vaf_max: float | None
    alt_reads_median: float | None
    alt_reads_min: int | None


class PlantingSummary(Contract):
    runs: list[str]
    samples: int
    variants: int
    complete: Proportion = Field(description="variants with status ok")
    by_status: dict[str, int]
    by_zygosity: list[ZygositySummary]


class FidelityResult(Contract):
    """Persisted as `fidelity.json`."""

    schema_version: int
    name: str
    baseline: str
    datasets: dict[str, str] = Field(description="id -> label")
    pass_only: bool
    recall_gap: float | None
    planting: PlantingSummary | None
    read_qc: list[ReadQc]
    callset_qc: list[CallsetQc]
    calling: list[CallingMetric]
    gaps: list[RecallGap]
    missing: list[str] = Field(description="raw outputs not found when assembling")
