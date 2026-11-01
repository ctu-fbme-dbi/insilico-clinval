"""Step 3b: read simulation with mason_simulator, true alignments and the planting check."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import ValidationError

from insilico_clinval.provenance import ManifestBuilder, sha256_file
from insilico_clinval.schemas.sample import SampleSpec
from insilico_clinval.schemas.simulation import SimulationProfile


class SimulationError(RuntimeError):
    """Simulation inputs are missing or inconsistent, or the planting check failed."""


def load_profile(path: Path) -> tuple[SimulationProfile, str]:
    """Load and validate a simulation profile; returns it with the file's SHA-256."""
    if not path.is_file():
        raise SimulationError(f"simulation profile not found: {path}")
    try:
        profile = SimulationProfile.model_validate(yaml.safe_load(path.read_text()))
    except (yaml.YAMLError, ValidationError) as exc:
        raise SimulationError(f"{path}: invalid simulation profile:\n{exc}") from exc
    return profile, sha256_file(path)


def load_truth_outputs(
    out_dir: Path, sample_id: str, manifest: ManifestBuilder
) -> tuple[SampleSpec, Path]:
    """SampleSpec and truth VCF written by build-truth; the VCF hash must match the manifest."""
    spec_path = out_dir / f"{sample_id}.sample.json"
    if not spec_path.is_file():
        raise SimulationError(f"{spec_path} not found (run `clinval build-truth` first)")
    spec = SampleSpec.model_validate_json(spec_path.read_text())
    truth = out_dir / spec.truth_vcf
    recorded = {o.role: o.sha256 for o in manifest.manifest.outputs}.get(f"truth_vcf:{sample_id}")
    if not truth.is_file() or recorded is None:
        raise SimulationError(f"truth VCF of {sample_id} missing or not in the run manifest")
    if sha256_file(truth) != recorded:
        raise SimulationError(f"{truth.name} differs from the file recorded in the run manifest")
    return spec, truth
