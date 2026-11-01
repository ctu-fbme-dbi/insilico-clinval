"""Enumerations from docs/schemas.md. All values are lower-case strings."""

from enum import Enum


class VariantCategory(str, Enum):
    STANDARD = "standard"
    CURATED = "curated"
    DIFFICULT = "difficult"
    SYNTHETIC_NOVEL = "synthetic_novel"


class VariantType(str, Enum):
    SNV = "snv"
    MNV = "mnv"
    INSERTION = "insertion"
    DELETION = "deletion"
    DELINS = "delins"


class Zygosity(str, Enum):
    HET = "het"
    HOM = "hom"
    HEMI = "hemi"


class Haplotype(str, Enum):
    ONE = "1"
    TWO = "2"
    BOTH = "both"


class Classification(str, Enum):
    PATHOGENIC = "pathogenic"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    UNCERTAIN_SIGNIFICANCE = "uncertain_significance"
    LIKELY_BENIGN = "likely_benign"
    BENIGN = "benign"


class SampleRole(str, Enum):
    POSITIVE = "positive"
    NEGATIVE = "negative"
    MIXED = "mixed"


class Sex(str, Enum):
    FEMALE = "female"
    MALE = "male"


class EvaluationLevel(str, Enum):
    CALLING = "calling"
    INTERPRETATION = "interpretation"
    BOTH = "both"


class Phase(str, Enum):
    """Relative phase of the members of a compound-heterozygous phase set."""

    CIS = "cis"
    TRANS = "trans"


class OutcomeClass(str, Enum):
    """Finding-level outcome, ordered along the pipeline (first failing stage wins)."""

    NOT_SIMULATED = "not_simulated"
    NOT_COVERED = "not_covered"
    MISMAPPED = "mismapped"
    NOT_CALLED = "not_called"
    FILTERED = "filtered"
    GENOTYPE_MISMATCH = "genotype_mismatch"
    NOT_ANNOTATED = "not_annotated"
    MISCLASSIFIED = "misclassified"
    REPORT_MISMATCH = "report_mismatch"
    CORRECT = "correct"
    FALSE_REPORT = "false_report"


class PlantingStatus(str, Enum):
    """Planting-completeness status of one planted variant in the simulated reads."""

    OK = "ok"
    MISSING_ON_HAPLOTYPE = "missing_on_haplotype"
    NO_READS = "no_reads"


class SimulationMode(str, Enum):
    REGION = "region"
    GENOME = "genome"
