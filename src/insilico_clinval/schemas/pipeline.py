"""PipelineConfig: the system under test (secondary + tertiary analysis), see configs/adapters/.

Separate from the suite: the same simulated data can be run through several pipelines (e.g. two
versions for a regression test). Paths are strings, resolved relative to the pipeline file.
"""

from typing import Any, Literal

from pydantic import Field, model_validator

from insilico_clinval.schemas.base import Contract
from insilico_clinval.schemas.enums import Classification
from insilico_clinval.schemas.runtime import RunnerName

# -- generic black-box command --------------------------------------------------------------


class CommandSpec(Contract):
    """An external command run as a black box (own runner, optionally its own image)."""

    tool_name: str = Field(pattern=r"^[A-Za-z0-9_.-]+$", description="name in the manifest")
    version: str = Field(description="recorded in the manifest (black boxes have no --version)")
    runner: RunnerName = "local"
    image: str | None = None
    argv: list[str] = Field(min_length=1, description="argv template with {placeholders}")

    @model_validator(mode="after")
    def _image(self) -> "CommandSpec":
        if self.runner == "docker" and self.image is None:
            raise ValueError(f"command {self.tool_name}: the docker runner needs an image")
        return self


# -- secondary ------------------------------------------------------------------------------


class SarekReference(Contract):
    """Either an iGenomes key or a custom FASTA (then igenomes_ignore is set)."""

    genome: str | None = Field(default=None, description="iGenomes key, e.g. GATK.GRCh38")
    fasta: str | None = None
    fasta_fai: str | None = None
    fasta_dict: str | None = Field(default=None, description="sarek parameter `dict`")
    bwa: str | None = Field(default=None, description="directory with a prebuilt BWA index")
    dbsnp: str | None = None
    dbsnp_tbi: str | None = None
    known_indels: str | None = None
    known_indels_tbi: str | None = None

    @model_validator(mode="after")
    def _one(self) -> "SarekReference":
        if (self.genome is None) == (self.fasta is None):
            raise ValueError("sarek reference: give exactly one of genome (iGenomes) or fasta")
        return self


class DeepVariantHaploid(Contract):
    """Haploid calling of sex chromosomes in XY samples (DeepVariant >= 1.6)."""

    contigs: list[str] = Field(default_factory=lambda: ["chrX", "chrY"], min_length=1)
    par_bed: str = Field(description="PAR regions, called as diploid (BED, e.g. GRCh38 PAR1/PAR2)")


class SarekConfig(Contract):
    release: str = Field(description="nf-core/sarek release tag, e.g. 3.10.0")
    profile: str = "docker"
    tools: list[str] = Field(default_factory=lambda: ["haplotypecaller"])
    skip_tools: list[str] = Field(default_factory=list)
    reference: SarekReference
    intervals: str | None = None
    wes: bool = False
    save_reference: bool = False
    variant_caller: Literal["haplotypecaller", "deepvariant"] = Field(
        default="haplotypecaller", description="caller whose VCF is collected"
    )
    params: dict[str, Any] = Field(
        default_factory=dict, description="further sarek parameters (passed in the params file)"
    )
    resource_limits: dict[str, Any] | None = Field(
        default=None, description="Nextflow process.resourceLimits, e.g. {cpus: 8, memory: 16.GB}"
    )
    extra_args: list[str] = Field(default_factory=list, description="extra nextflow arguments")
    deepvariant_haploid: DeepVariantHaploid | None = Field(
        default=None,
        description="call these contigs as haploid in XY samples (sarek passes no such option)",
    )

    @model_validator(mode="after")
    def _haploid_needs_deepvariant(self) -> "SarekConfig":
        if self.deepvariant_haploid is not None and "deepvariant" not in self.tools:
            raise ValueError("deepvariant_haploid needs `deepvariant` in sarek.tools")
        return self


class SecondaryOutputs(Contract):
    """Output path templates (relative to the adapter's output directory; {sample})."""

    alignment: str | None = None
    raw_vcf: str
    filtered_vcf: str | None = None


class CommandSecondary(Contract):
    command: CommandSpec
    outputs: SecondaryOutputs


class PrecomputedSample(Contract):
    alignment: str | None = None
    raw_vcf: str
    filtered_vcf: str | None = None


class PrecomputedSecondary(Contract):
    tool_name: str = "precomputed"
    version: str
    samples: dict[str, PrecomputedSample]


class SecondaryConfig(Contract):
    adapter: Literal["sarek", "command", "precomputed"]
    sarek: SarekConfig | None = None
    command: CommandSecondary | None = None
    precomputed: PrecomputedSecondary | None = None

    @model_validator(mode="after")
    def _section(self) -> "SecondaryConfig":
        _require_section(self, self.adapter, ("sarek", "command", "precomputed"))
        return self


# -- tertiary -------------------------------------------------------------------------------


class NirvanaData(Contract):
    cache: str = Field(description="Nirvana cache prefix (e.g. …/GRCh38/Cache/Both)")
    supplementary: str = Field(description="SupplementaryAnnotation directory for GRCh38")
    reference: str = Field(description="Homo_sapiens.GRCh38.Nirvana.dat")


class VepData(Contract):
    """Ensembl VEP inputs for BIAS-2015's VEP mode (BIAS doc/vep_setup.md, hg38)."""

    cache_dir: str = Field(description="VEP cache root (contains homo_sapiens/<version>_GRCh38)")
    cache_version: str = Field(description="VEP cache release, e.g. 115 (must match the image)")
    fasta: str = Field(description="VEP's GRCh38 FASTA (bgzipped, indexed)")
    gnomad_exomes: str = Field(description="gnomAD exomes sites VCF with AF, AF_joint, AC, AN")
    clinvar: str = Field(description="ClinVar GRCh38 VCF")
    clinvar_release: str = Field(description="ClinVar release of the VCF, e.g. 2026-02-01")
    alphamissense: str = Field(description="AlphaMissense hg38 TSV for the VEP plugin")
    revel: str = Field(description="REVEL GRCh38 TSV for the VEP plugin")
    fork: int = Field(default=4, ge=1, description="VEP --fork")


class Bias2015Config(Contract):
    version: str = Field(description="BIAS-2015 release, e.g. v3.0.0 (must match the image)")
    annotator: Literal["nirvana", "vep"] = Field(
        default="nirvana", description="annotator of BIAS's input JSON"
    )
    data_dir: str = Field(description="BIAS hg38 data bundle directory")
    data_release: str = Field(description="data bundle release, e.g. 2026.03.01")
    nirvana: NirvanaData | None = None
    vep: VepData | None = None
    skip_codes: list[str] = Field(
        default_factory=list, description="ACMG codes BIAS must not assign, e.g. [PP5, BP6]"
    )

    @model_validator(mode="after")
    def _annotator(self) -> "Bias2015Config":
        _require_section(self, self.annotator, ("nirvana", "vep"))
        return self


class ParserSpec(Contract):
    """How to read a black-box classifier output into CanonicalFindings."""

    format: Literal["tsv", "csv", "json"]
    records: str | None = Field(default=None, description="json: key of the record list")
    columns: dict[
        Literal[
            "chrom", "pos", "ref", "alt", "classification", "criteria", "gene", "hgvs_c",
            "hgvs_p", "reported",
        ],
        str,
    ]  # fmt: skip
    classification_map: dict[str, Classification]
    criteria_separator: str = ","
    reported_true: list[str] = Field(default_factory=lambda: ["true", "yes", "1", "reported"])

    @model_validator(mode="after")
    def _required(self) -> "ParserSpec":
        missing = [
            c for c in ("chrom", "pos", "ref", "alt", "classification") if c not in self.columns
        ]
        if missing:
            raise ValueError(f"parser.columns misses {missing}")
        return self


class CommandTertiary(Contract):
    command: CommandSpec
    output: str = Field(description="output file template relative to the work directory")
    parser: ParserSpec


class TertiaryConfig(Contract):
    adapter: Literal["bias2015", "command", "none"]
    input: Literal["filtered", "raw"] = "filtered"
    pass_only: bool = True
    regions_bed: str | None = Field(default=None, description="restrict to reportable regions")
    report_classifications: list[Classification] = Field(
        default_factory=lambda: [Classification.PATHOGENIC, Classification.LIKELY_PATHOGENIC]
    )
    bias2015: Bias2015Config | None = None
    command: CommandTertiary | None = None

    @model_validator(mode="after")
    def _section(self) -> "TertiaryConfig":
        _require_section(self, self.adapter, ("bias2015", "command"))
        return self


class PipelineConfig(Contract):
    schema_version: Literal[1]
    name: str = Field(pattern=r"^[A-Za-z0-9_.-]+$", description="directory name under sut/")
    description: str = ""
    secondary: SecondaryConfig
    tertiary: TertiaryConfig


def _require_section(model: Contract, adapter: str, sections: tuple[str, ...]) -> None:
    for name in sections:
        present = getattr(model, name) is not None
        if name == adapter and not present:
            raise ValueError(f"adapter {adapter} needs a `{name}` section")
        if name != adapter and present:
            raise ValueError(f"section `{name}` given but adapter is {adapter}")
