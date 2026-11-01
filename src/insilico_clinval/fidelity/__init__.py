"""
Fidelity of simulated vs real data: planting completeness across runs, read-level QC,
call-set statistics and stratified calling performance (hap.py) of a simulated and a real genome
processed by the same secondary pipeline.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import yaml
from pydantic import ValidationError

from insilico_clinval.schemas.fidelity import FidelityConfig, FidelityDataset


class FidelityConfigError(ValueError):
    """The fidelity configuration is invalid or references missing files."""


@dataclass(frozen=True)
class LoadedFidelity:
    config: FidelityConfig
    path: Path
    sha256: str

    def resolve(self, value: str) -> Path:
        p = Path(value)
        return (p if p.is_absolute() else self.path.parent / p).resolve()

    @property
    def reference(self) -> Path:
        return self.resolve(self.config.reference)

    @property
    def confident_bed(self) -> Path:
        return self.resolve(self.config.confident_bed)

    @property
    def qc_regions(self) -> Path:
        cfg = self.config
        return self.resolve(cfg.qc_regions if cfg.qc_regions is not None else cfg.confident_bed)

    def strata(self) -> list[tuple[str, Path]]:
        return [(name, self.resolve(p)) for name, p in sorted(self.config.stratifications.items())]

    def dataset(self, dataset_id: str) -> FidelityDataset:
        for d in self.config.datasets:
            if d.id == dataset_id:
                return d
        raise KeyError(dataset_id)


def load_fidelity(path: Path) -> LoadedFidelity:
    if not path.is_file():
        raise FidelityConfigError(f"fidelity config not found: {path}")
    raw = path.read_bytes()
    try:
        config = FidelityConfig.model_validate(yaml.safe_load(raw))
    except (yaml.YAMLError, ValidationError) as exc:
        raise FidelityConfigError(f"{path}: invalid fidelity config:\n{exc}") from exc
    return LoadedFidelity(
        config=config, path=path.resolve(), sha256=hashlib.sha256(raw).hexdigest()
    )


def referenced_files(fid: LoadedFidelity) -> list[tuple[str, Path]]:
    """Input files of the configuration as ``(role, path)``, in configuration order."""
    cfg = fid.config
    items = [("reference", fid.reference), ("confident_bed", fid.confident_bed)]
    if cfg.qc_regions is not None:
        items.append(("qc_regions", fid.qc_regions))
    items += [(f"stratum:{name}", path) for name, path in fid.strata()]
    for d in cfg.datasets:
        items += [
            (f"alignment:{d.id}", fid.resolve(d.alignment)),
            (f"vcf:{d.id}", fid.resolve(d.vcf)),
            (f"truth_vcf:{d.id}", fid.resolve(d.truth_vcf)),
        ]
    if cfg.happy.vcfeval_template is not None:
        items.append(("vcfeval_template", fid.resolve(cfg.happy.vcfeval_template)))
    items += [(f"planting_run:{i}", fid.resolve(p)) for i, p in enumerate(cfg.planting)]
    return items


def missing_files(fid: LoadedFidelity) -> list[str]:
    return [f"{role}: {path}" for role, path in referenced_files(fid) if not path.exists()]
