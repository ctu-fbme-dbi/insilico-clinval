"""Steps 4-5: adapters that run the system under test as a black box and normalise its outputs.

Secondary analysis (FASTQ -> alignment + VCF): ``sarek``, ``command``, ``precomputed``.
Tertiary analysis (VCF -> classified, reported findings): ``bias2015``, ``command``, ``none``.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import yaml
from pydantic import ValidationError

from insilico_clinval.schemas.enums import Sex
from insilico_clinval.schemas.pipeline import PipelineConfig


class PipelineConfigError(ValueError):
    """The pipeline config is invalid or references missing files."""


class SutFailure(RuntimeError):
    """The system under test failed or did not produce its declared outputs."""


@dataclass(frozen=True)
class SampleInput:
    """Simulated reads of one sample, the input of secondary analysis."""

    sample_id: str
    sex: Sex
    fastq_r1: Path
    fastq_r2: Path


@dataclass(frozen=True)
class LoadedPipeline:
    config: PipelineConfig
    path: Path
    sha256: str

    def resolve(self, value: str) -> Path:
        p = Path(value)
        return (p if p.is_absolute() else self.path.parent / p).resolve()


def load_pipeline(path: Path) -> LoadedPipeline:
    if not path.is_file():
        raise PipelineConfigError(f"pipeline config not found: {path}")
    raw = path.read_bytes()
    try:
        config = PipelineConfig.model_validate(yaml.safe_load(raw))
    except (yaml.YAMLError, ValidationError) as exc:
        raise PipelineConfigError(f"{path}: invalid pipeline config:\n{exc}") from exc
    return LoadedPipeline(
        config=config, path=path.resolve(), sha256=hashlib.sha256(raw).hexdigest()
    )


def referenced_files(pipeline: LoadedPipeline) -> list[tuple[str, Path]]:
    """Files and directories the pipeline config refers to, as ``(role, path)``."""
    cfg = pipeline.config
    items: list[tuple[str, Path]] = []
    sec = cfg.secondary
    if sec.sarek is not None:
        ref = sec.sarek.reference
        for key in (
            "fasta", "fasta_fai", "fasta_dict", "bwa", "dbsnp", "dbsnp_tbi", "known_indels",
            "known_indels_tbi",
        ):  # fmt: skip
            value = getattr(ref, key)
            if value is not None:
                items.append((f"sarek_{key}", pipeline.resolve(value)))
        if sec.sarek.intervals is not None:
            items.append(("sarek_intervals", pipeline.resolve(sec.sarek.intervals)))
    if sec.precomputed is not None:
        for sid, files in sorted(sec.precomputed.samples.items()):
            for key in ("alignment", "raw_vcf", "filtered_vcf"):
                value = getattr(files, key)
                if value is not None:
                    items.append((f"precomputed_{key}:{sid}", pipeline.resolve(value)))
    ter = cfg.tertiary
    if ter.regions_bed is not None:
        items.append(("tertiary_regions_bed", pipeline.resolve(ter.regions_bed)))
    if ter.bias2015 is not None:
        b = ter.bias2015
        items.append(("bias2015_data_dir", pipeline.resolve(b.data_dir)))
        if b.nirvana is not None:
            items += [
                ("nirvana_cache", pipeline.resolve(b.nirvana.cache).parent),
                ("nirvana_supplementary", pipeline.resolve(b.nirvana.supplementary)),
                ("nirvana_reference", pipeline.resolve(b.nirvana.reference)),
            ]
        if b.vep is not None:
            v = b.vep
            items.append(("vep_cache_dir", pipeline.resolve(v.cache_dir)))
            for key in ("fasta", "gnomad_exomes", "clinvar", "alphamissense", "revel"):
                items.append((f"vep_{key}", pipeline.resolve(getattr(v, key))))
    return items


def missing_files(pipeline: LoadedPipeline) -> list[str]:
    return [f"{role}: {path}" for role, path in referenced_files(pipeline) if not path.exists()]
