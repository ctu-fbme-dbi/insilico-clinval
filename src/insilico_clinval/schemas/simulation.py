"""Simulation contracts: profile (YAML input), simulation spec and planting check (outputs)."""

from typing import Literal

from pydantic import Field, model_validator

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import PlantingStatus, SimulationMode, Zygosity

# Mason's minimum --chunk-size (mason_options.cpp, seqan-v2.4.0).
MASON_MIN_CHUNK_SIZE = 65536


class FragmentModel(Contract):
    model: Literal["normal", "uniform"] = "normal"
    mean: int = Field(ge=1, description="mean fragment (insert) size, bp")
    sd: int = Field(ge=0, description="standard deviation (normal model)")
    min: int | None = Field(default=None, ge=1, description="minimum size (uniform model)")
    max: int | None = Field(default=None, ge=1, description="maximum size (uniform model)")

    @model_validator(mode="after")
    def _uniform(self) -> "FragmentModel":
        if self.model == "uniform" and (self.min is None or self.max is None):
            raise ValueError("the uniform fragment model needs min and max")
        if self.min is not None and self.max is not None and self.max < self.min:
            raise ValueError("fragment max < min")
        return self

    @property
    def longest(self) -> int:
        """Practical upper bound of the fragment size (max, else mean + 4 sd)."""
        return self.max if self.max is not None else self.mean + 4 * self.sd


class IlluminaErrors(Contract):
    """Optional Mason Illumina error-model parameters; null = Mason's default."""

    prob_mismatch: float | None = Field(default=None, ge=0, le=1)
    prob_mismatch_scale: float | None = Field(default=None, ge=0)
    prob_mismatch_begin: float | None = Field(default=None, ge=0, le=1)
    prob_mismatch_end: float | None = Field(default=None, ge=0, le=1)
    prob_insert: float | None = Field(default=None, ge=0, le=1)
    prob_deletion: float | None = Field(default=None, ge=0, le=1)


class SimulationProfile(Contract):
    """Read-simulation parameters (configs/profiles/*.yaml)."""

    schema_version: Literal[1]
    name: str
    technology: Literal["illumina"]
    read_length: int = Field(ge=20)
    fragment: FragmentModel
    coverage: float = Field(gt=0, description="mean depth of a diploid locus (both haplotypes)")
    threads: int = Field(
        default=1, ge=1, description="Mason threads; part of the determinism key (one RNG each)"
    )
    chunk_size: int = Field(default=MASON_MIN_CHUNK_SIZE, ge=MASON_MIN_CHUNK_SIZE)
    illumina: IlluminaErrors = Field(default_factory=IlluminaErrors)

    @model_validator(mode="after")
    def _fragment(self) -> "SimulationProfile":
        if self.fragment.mean < self.read_length:
            raise ValueError("fragment mean must be at least the read length")
        return self


class SegmentRecord(Contract):
    """One simulated reference segment (1-based closed genome interval)."""

    name: str
    chrom: str
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    ploidy: Literal[1, 2]


class MasonRun(Contract):
    """One mason_simulator invocation (one ploidy group, or one batch of it, of one sample)."""

    group: Literal["diploid", "haploid"]
    name: str | None = Field(
        default=None,
        description="batch name (`diploid-1`, ...) with simulation.mason_batch_bp; null otherwise",
    )
    segment_names: list[str] | None = Field(
        default=None, description="segments of the batch; null without simulation.mason_batch_bp"
    )
    ploidy: Literal[1, 2]
    segments: int
    length: int = Field(description="total segment length, bp")
    fragments: int
    seed: int
    vcf: str | None = Field(description="Mason input VCF, relative to the run directory")
    vcf_records: int
    args: list[str] = Field(description="mason_simulator arguments (paths relative to its cwd)")


class SimulationSpec(Contract):
    """Persisted as `<sample_id>.simulation.json`."""

    schema_version: int
    sample_id: str
    mode: SimulationMode
    profile: SimulationProfile
    profile_sha256: str
    flank_bp: int | None
    segments: str = Field(description="segments.tsv, relative to the run directory")
    runs: list[MasonRun]
    decomposed: list[str] = Field(description="planted MNV/delins ids split for Mason")
    background_dropped: int
    background_dropped_log: str | None
    fastq_r1: str
    fastq_r2: str
    true_bam: str
    planting_check: str | None = None
    planting_complete: bool | None = None


class HaplotypeReads(Contract):
    haplotype: int = Field(description="Mason haplotype (oH tag), 1-based")
    reads: int
    alt: int
    ref: int
    other: int


class VariantPlanting(Contract):
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    zygosity: Zygosity
    carrying_haplotypes: list[int]
    expected_vaf: float
    spanning_reads: int
    alt_reads: int
    ref_reads: int
    other_reads: int
    vaf: float | None
    by_haplotype: list[HaplotypeReads]
    status: PlantingStatus


class PlantingCheck(Contract):
    """Persisted as `<sample_id>.planting_check.json`."""

    schema_version: int
    sample_id: str
    min_alt_reads: int
    window_pad: int
    complete: bool
    variants: list[VariantPlanting]
