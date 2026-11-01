"""Runner abstraction: every external tool is executed ``local``ly or in ``docker``.

The runner only builds the command line; execution, return-code checks, stderr capture and
resource measurement are done by the single subprocess wrapper ``_exec.run``. The
:class:`Toolbox` resolves the runner of a tool from a :class:`RuntimeConfig`, captures the tool
version through the same runner and records version, container digest and step in the manifest.

Docker runs use ``docker run --rm --user UID:GID``; every mount is bound at its resolved path, so
paths in the argument list are valid inside the container unchanged. Note: for the docker runner,
the peak RSS measured by ``_exec.run`` is that of the docker client, not of the tool, and is
therefore not recorded.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import yaml
from pydantic import ValidationError

from insilico_clinval._exec import ExecResult, ToolError, run
from insilico_clinval.schemas.runtime import KNOWN_TOOLS, RuntimeConfig

if TYPE_CHECKING:
    from insilico_clinval.provenance import ManifestBuilder

log = logging.getLogger(__name__)


class RuntimeConfigError(ValueError):
    """The runtime configuration file is invalid."""


def load_runtime(path: Path | None) -> RuntimeConfig:
    """Load a runtime YAML; ``None`` gives the default (every tool local)."""
    if path is None:
        return RuntimeConfig()
    if not path.is_file():
        raise RuntimeConfigError(f"runtime config not found: {path}")
    try:
        data = yaml.safe_load(path.read_text())
        return RuntimeConfig.model_validate(data)
    except (yaml.YAMLError, ValidationError) as exc:
        raise RuntimeConfigError(f"{path}: invalid runtime config:\n{exc}") from exc


def collapse_mounts(paths: Sequence[Path]) -> list[Path]:
    """Resolved, de-duplicated mount points; paths inside another mount are dropped."""
    resolved = sorted({p.resolve() for p in paths}, key=lambda p: (len(p.parts), str(p)))
    out: list[Path] = []
    for p in resolved:
        if not any(p == m or m in p.parents for m in out):
            out.append(p)
    return sorted(out, key=str)


# Version commands that differ from ``<tool> --version``.
DEFAULT_VERSION_ARGS: dict[str, list[str]] = {"nextflow": ["-version"], "vep": ["--help"]}
# Tools whose version output starts with a banner: take the first line mentioning the marker.
BANNER_TOOLS: dict[str, str] = {"nextflow": "version", "vep": "ensembl-vep"}


@dataclass(frozen=True)
class BlackBox:
    """A command of the system under test, run as a black box (argv given in full)."""

    runner: str
    image: str | None
    version: str


def version_line(text: str, *, banner: str | None = None) -> str:
    """First non-empty line; with ``banner``, the first line containing that marker."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if banner:
        for line in lines:
            if banner in line.lower():
                return line
    return lines[0] if lines else "unknown"


class Toolbox:
    """Executes external tools according to a runtime configuration."""

    def __init__(self, runtime: RuntimeConfig | None = None):
        self.runtime = runtime if runtime is not None else RuntimeConfig()
        self._versions: dict[str, str] = {}
        self._digests: dict[str, str | None] = {}
        self._black_boxes: dict[str, BlackBox] = {}

    def register_black_box(
        self, name: str, *, runner: str, image: str | None, version: str
    ) -> None:
        """Register a black-box command (its version is declared, not queried)."""
        if name in KNOWN_TOOLS:
            raise ToolError(f"black-box name {name!r} clashes with a framework tool")
        self._black_boxes[name] = BlackBox(runner=runner, image=image, version=version)
        self._versions[name] = version

    # -- command construction -------------------------------------------------------------

    def runner(self, tool: str) -> str:
        if tool in self._black_boxes:
            return self._black_boxes[tool].runner
        return self.runtime.runner_of(tool)

    def executable(self, tool: str) -> str:
        spec = self.runtime.tools.get(tool)
        return spec.executable if spec is not None and spec.executable else tool

    def image(self, tool: str) -> str:
        if tool in self._black_boxes:
            image = self._black_boxes[tool].image
        else:
            spec = self.runtime.tools.get(tool)
            image = spec.image if spec is not None else None
        if image is None:
            raise ToolError(f"tool {tool}: the docker runner needs an image")
        return image

    def command(
        self,
        tool: str,
        args: Sequence[str],
        *,
        mounts: Sequence[Path] = (),
        cwd: Path | None = None,
        env: Mapping[str, str] | None = None,
    ) -> list[str]:
        """The argv that runs ``tool args`` with the configured runner.

        For a black box, ``args`` is the complete command (executable included). ``env`` becomes
        ``-e`` options for docker; a local tool gets it from ``run`` (not part of the argv).
        """
        if tool in self._black_boxes:
            head: list[str] = []
        elif tool in KNOWN_TOOLS:
            head = [self.executable(tool)]
        else:
            raise ToolError(f"unknown tool {tool!r} (known: {', '.join(KNOWN_TOOLS)})")
        if self.runner(tool) == "local":
            return [*head, *args]
        docker = self.runtime.docker
        argv = [docker.binary, "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}"]
        argv += list(docker.extra_args)
        all_mounts = [*mounts, cwd] if cwd is not None else list(mounts)
        for m in collapse_mounts(all_mounts):
            argv += ["-v", f"{m}:{m}"]
        if cwd is not None:
            argv += ["-w", str(cwd.resolve())]
        for key, value in sorted((env or {}).items()):
            argv += ["-e", f"{key}={value}"]
        return [*argv, self.image(tool), *head, *args]

    # -- provenance -----------------------------------------------------------------------

    def version(self, tool: str) -> str:
        """Version line of ``<tool> --version`` or the configured version command (cached)."""
        if tool not in self._versions:
            spec = self.runtime.tools.get(tool)
            args = spec.version_args if spec is not None and spec.version_args else None
            args = args or DEFAULT_VERSION_ARGS.get(tool, ["--version"])
            result = run(self.command(tool, args))
            text = result.stdout if result.stdout.strip() else result.stderr
            self._versions[tool] = version_line(text, banner=BANNER_TOOLS.get(tool))
        return self._versions[tool]

    def container_digest(self, tool: str) -> str | None:
        """``sha256:…`` of the image (pinned digest, else the local image id); None if local."""
        if self.runner(tool) == "local":
            return None
        if tool not in self._digests:
            image = self.image(tool)
            if "@sha256:" in image:
                self._digests[tool] = image.split("@", 1)[1]
            else:
                result = run(
                    [self.runtime.docker.binary, "image", "inspect", "--format", "{{.Id}}", image]
                )
                self._digests[tool] = result.stdout.strip() or None
        return self._digests[tool]

    def record(self, tool: str, manifest: ManifestBuilder) -> None:
        manifest.add_tool(tool, self.version(tool), self.container_digest(tool))
        if self.runner(tool) == "docker" and "@sha256:" not in self.image(tool):
            manifest.warn(
                f"tool {tool}: docker image {self.image(tool)} is not pinned by digest "
                f"(image id {self.container_digest(tool)})"
            )

    # -- execution ------------------------------------------------------------------------

    def run(
        self,
        tool: str,
        args: Sequence[str],
        *,
        step: str,
        manifest: ManifestBuilder | None = None,
        mounts: Sequence[Path] = (),
        cwd: Path | None = None,
        stderr_path: Path | None = None,
        env: Mapping[str, str] | None = None,
    ) -> ExecResult:
        """Run a tool, recording version, container digest and the step in ``manifest``.

        ``stderr_path``: write the tool's stderr there while it runs (kept on failure).
        ``env``: extra environment variables for the tool (local or inside the container).
        """
        argv = self.command(tool, args, mounts=mounts, cwd=cwd, env=env)
        local_env = env if self.runner(tool) == "local" else None
        if manifest is not None:
            self.record(tool, manifest)
        result = run(argv, cwd=cwd, check=False, stderr_path=stderr_path, env=local_env)
        if manifest is not None:
            manifest.add_step(
                name=step,
                started_at=result.started_at,
                wall_s=result.wall_s,
                peak_rss_mb=result.peak_rss_mb if self.runner(tool) == "local" else None,
                exit_code=result.returncode,
            )
        if result.returncode != 0:
            raise ToolError(
                f"{tool} exited with status {result.returncode}: {' '.join(argv)}\n"
                f"{result.stderr.strip()[-4000:]}"
            )
        return result
