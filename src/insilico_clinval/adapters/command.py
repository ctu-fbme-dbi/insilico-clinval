"""Generic black-box command: placeholder substitution and execution through the runner."""

from __future__ import annotations

import string
from collections.abc import Mapping, Sequence
from pathlib import Path

from insilico_clinval._runner import Toolbox
from insilico_clinval.adapters import PipelineConfigError
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas.pipeline import CommandSpec

_FORMATTER = string.Formatter()


def substitute(template: str, values: Mapping[str, str]) -> str:
    """Replace ``{name}`` placeholders; unknown placeholders are configuration errors."""
    names = [field for _, field, _, _ in _FORMATTER.parse(template) if field is not None]
    unknown = sorted(set(names) - set(values))
    if unknown:
        raise PipelineConfigError(
            f"unknown placeholder(s) {unknown} in {template!r}; available: {sorted(values)}"
        )
    return template.format_map(dict(values))


def run_black_box(
    toolbox: Toolbox,
    spec: CommandSpec,
    values: Mapping[str, str],
    *,
    workdir: Path,
    step: str,
    manifest: ManifestBuilder | None,
    mounts: Sequence[Path] = (),
) -> None:
    """Run ``spec.argv`` (placeholders substituted) with its own runner, in ``workdir``."""
    workdir.mkdir(parents=True, exist_ok=True)
    toolbox.register_black_box(
        spec.tool_name, runner=spec.runner, image=spec.image, version=spec.version
    )
    argv = [substitute(arg, values) for arg in spec.argv]
    toolbox.run(
        spec.tool_name, argv, step=step, manifest=manifest, mounts=[*mounts, workdir], cwd=workdir
    )
