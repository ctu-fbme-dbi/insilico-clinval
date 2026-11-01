"""Loading and validation of the YAML suite configuration.

Relative paths in the suite are resolved against the directory containing the suite file.
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from insilico_clinval.schemas.suite import ClinGenRule, ClinVarRule, GiabRule, SuiteConfig


class SuiteError(ValueError):
    """The suite file is invalid or references missing inputs."""


class _UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that rejects duplicate mapping keys instead of silently overwriting."""


def _construct_mapping(loader: _UniqueKeyLoader, node: yaml.MappingNode) -> dict[Any, Any]:
    mapping: dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node)
        if key in mapping:
            raise SuiteError(f"duplicate key {key!r} at line {key_node.start_mark.line + 1}")
        mapping[key] = loader.construct_object(value_node, deep=True)
    return mapping


_UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_mapping,
)


@dataclass(frozen=True)
class LoadedSuite:
    config: SuiteConfig
    path: Path
    sha256: str

    @property
    def base_dir(self) -> Path:
        return self.path.parent

    def resolve(self, value: str) -> Path:
        p = Path(value)
        return p if p.is_absolute() else (self.base_dir / p)

    def output_dir(self, override: Path | None) -> Path:
        """Run directory with ``..`` collapsed (lexically, so a ``runs`` symlink is kept)."""
        path = override if override is not None else self.resolve(self.config.output.directory)
        return Path(os.path.normpath(path))


def load_suite(path: Path) -> LoadedSuite:
    """Parse and validate a suite file (schema and cross-references, not file existence)."""
    if not path.is_file():
        raise SuiteError(f"suite file not found: {path}")
    raw = path.read_bytes()
    try:
        data = yaml.load(raw, Loader=_UniqueKeyLoader)  # safe: SafeLoader subclass
    except yaml.YAMLError as exc:
        raise SuiteError(f"{path}: invalid YAML: {exc}") from exc
    try:
        config = SuiteConfig.model_validate(data)
    except ValidationError as exc:
        raise SuiteError(f"{path}: invalid suite:\n{exc}") from exc
    suite = LoadedSuite(config=config, path=path.resolve(), sha256=hashlib.sha256(raw).hexdigest())
    _check_background_conflicts(suite)
    return suite


def _check_background_conflicts(suite: LoadedSuite) -> None:
    """GIAB rules must not plant background variants into samples that carry that background."""
    cfg = suite.config
    if cfg.background is None or cfg.sources.giab is None:
        return
    same = (
        suite.resolve(cfg.background.vcf).resolve()
        == suite.resolve(cfg.sources.giab.truth_vcf).resolve()
    )
    if not same:
        return
    with_bg = {s.sample_id for s in cfg.samples if cfg.uses_background(s)}
    errors = [
        f"rule {r.id}: selects from the background VCF ({cfg.background.name}) into samples "
        f"{sorted(with_bg.intersection(r.samples))} that already carry it as background; the "
        "background itself is the calling truth of those samples"
        for r in cfg.selection
        if isinstance(r, GiabRule) and with_bg.intersection(r.samples)
    ]
    if errors:
        raise SuiteError(f"{suite.path}: invalid suite:\n" + "\n".join(errors))


def required_inputs(suite: LoadedSuite) -> list[tuple[str, Path]]:
    """All input files the suite needs, as ``(role, path)`` in a stable order."""
    cfg = suite.config
    fasta = suite.resolve(cfg.reference.fasta)
    items: list[tuple[str, Path]] = [
        ("reference_fasta", fasta),
        ("reference_fai", Path(f"{fasta}.fai")),
    ]
    uses_erepo = any(isinstance(r, ClinGenRule) for r in cfg.selection)
    uses_clinvar = any(isinstance(r, ClinVarRule) for r in cfg.selection)
    uses_giab = any(isinstance(r, GiabRule) for r in cfg.selection)
    erepo = cfg.sources.clingen_erepo
    if erepo is not None and uses_erepo:
        items.append(("clingen_erepo", suite.resolve(erepo.file)))
        if erepo.held_out_reference is not None:
            items.append(("clingen_erepo_held_out", suite.resolve(erepo.held_out_reference.file)))
    if cfg.sources.clinvar is not None and uses_clinvar:
        items.append(("clinvar_vcf", suite.resolve(cfg.sources.clinvar.vcf)))
    if cfg.sources.giab is not None and uses_giab:
        items.append(("giab_truth_vcf", suite.resolve(cfg.sources.giab.truth_vcf)))
        items.append(("giab_confident_bed", suite.resolve(cfg.sources.giab.confident_bed)))
    if cfg.background is not None and any(cfg.uses_background(s) for s in cfg.samples):
        bg_vcf = suite.resolve(cfg.background.vcf)
        items.append(("background_vcf", bg_vcf))
        items.append(("background_vcf_index", Path(f"{bg_vcf}.tbi")))
        items.append(("background_bed", suite.resolve(cfg.background.bed)))
    for name in sorted(cfg.sources.strata):
        items.append((f"stratum:{name}", suite.resolve(cfg.sources.strata[name])))
    for rule in cfg.selection:
        if isinstance(rule, (ClinGenRule, ClinVarRule, GiabRule)) and rule.region_bed is not None:
            items.append((f"region_bed:{rule.id}", suite.resolve(rule.region_bed)))
    return items


def check_inputs_exist(suite: LoadedSuite) -> None:
    missing = [f"{role}: {path}" for role, path in required_inputs(suite) if not path.is_file()]
    if missing:
        raise SuiteError("missing required inputs:\n  " + "\n  ".join(missing))
    expected = suite.config.reference.sha256
    if expected is not None:
        from insilico_clinval.provenance import sha256_file

        fasta = suite.resolve(suite.config.reference.fasta)
        actual = sha256_file(fasta)
        if actual != expected:
            raise SuiteError(f"reference FASTA SHA-256 mismatch: expected {expected}, got {actual}")
