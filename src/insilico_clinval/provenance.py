"""Run manifest create/extend/write, file hashing and framework provenance."""

from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from insilico_clinval import __version__
from insilico_clinval._exec import ToolError, run
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.manifest import (
    ConfigRef,
    DatabaseRef,
    FileRef,
    FrameworkInfo,
    RunManifest,
    StepRecord,
    ToolRef,
)

MANIFEST_NAME = "run_manifest.json"
_CHUNK = 1 << 20


class ManifestError(RuntimeError):
    """The run directory holds a manifest that belongs to a different configuration."""


def sha256_file(path: Path) -> str:
    """SHA-256 of a file, streamed in 1 MiB chunks."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(_CHUNK):
            h.update(chunk)
    return h.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# Libraries whose version changes output bytes: pysam (htslib: BAM, VCF, BGZF), isal (FASTQ gzip).
OUTPUT_PACKAGES = ("isal", "pysam")


def python_packages() -> dict[str, str]:
    out: dict[str, str] = {}
    for name in OUTPUT_PACKAGES:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = "not installed"
    return out


def framework_info() -> FrameworkInfo:
    """Framework version, git commit and dirty flag of the checkout, output-shaping libraries."""
    src_dir = Path(__file__).resolve().parent
    try:
        commit = run(["git", "rev-parse", "HEAD"], cwd=src_dir).stdout.strip() or None
        status = run(["git", "status", "--porcelain"], cwd=src_dir).stdout
        dirty: bool | None = bool(status.strip())
    except ToolError:
        commit, dirty = None, None
    return FrameworkInfo(
        version=__version__, git_commit=commit, dirty=dirty, python_packages=python_packages()
    )


def _self_peak_rss_mb() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(ru / (1024 * 1024) if sys.platform == "darwin" else ru / 1024, 1)


class ManifestBuilder:
    """Accumulates provenance for one command and merges it into the run directory manifest."""

    def __init__(self, run_dir: Path, *, config_path: Path, config_sha256: str, seed: int):
        self.run_dir = run_dir
        self.path = run_dir / MANIFEST_NAME
        if self.path.exists():
            existing = RunManifest.model_validate_json(self.path.read_text())
            if existing.config.sha256 != config_sha256 or existing.seed != seed:
                raise ManifestError(
                    f"{self.path} was written for a different suite config or seed; "
                    "use a fresh output directory"
                )
            self.manifest = existing
            self._merge_python_packages()
        else:
            created = _now()
            self.manifest = RunManifest(
                schema_version=SCHEMA_VERSION,
                run_id=f"{created.replace(':', '').replace('-', '')[:15]}-{config_sha256[:8]}",
                created_at=created,
                framework=framework_info(),
                config=ConfigRef(path=str(config_path), sha256=config_sha256),
                seed=seed,
                inputs=[],
                outputs=[],
                tools=[],
                databases=[],
                steps=[],
                warnings=[],
            )

    def _merge_python_packages(self) -> None:
        """Add this process's library versions to an existing manifest; warn when they differ."""
        recorded = self.manifest.framework.python_packages
        current = python_packages()
        for name, ver in current.items():
            if name in recorded and recorded[name] != ver:
                self.warn(
                    f"{name} {ver} differs from {recorded[name]} recorded earlier in this run"
                )
        merged = {**current, **recorded}
        if merged != recorded:
            self.manifest.framework = self.manifest.framework.model_copy(
                update={"python_packages": dict(sorted(merged.items()))}
            )

    @staticmethod
    def _upsert(items: list[FileRef], ref: FileRef) -> None:
        for i, existing in enumerate(items):
            if existing.role == ref.role and existing.path == ref.path:
                items[i] = ref
                return
        items.append(ref)

    def add_input(self, role: str, path: Path) -> None:
        ref = FileRef(role=role, path=str(path), sha256=sha256_file(path))
        self._upsert(self.manifest.inputs, ref)

    def add_output(self, role: str, path: Path) -> None:
        """Record an output file; the path is stored relative to the run directory."""
        rel = path.resolve().relative_to(self.run_dir.resolve())
        self._upsert(
            self.manifest.outputs, FileRef(role=role, path=str(rel), sha256=sha256_file(path))
        )

    def add_tool(self, name: str, version: str, container_digest: str | None = None) -> None:
        tools = self.manifest.tools
        tools[:] = [t for t in tools if t.name != name]
        tools.append(ToolRef(name=name, version=version, container_digest=container_digest))
        tools.sort(key=lambda t: t.name)

    def add_database(self, name: str, release: str, sha256: str | None) -> None:
        dbs = self.manifest.databases
        dbs[:] = [d for d in dbs if d.name != name]
        dbs.append(DatabaseRef(name=name, release=release, sha256=sha256))
        dbs.sort(key=lambda d: d.name)

    def add_step(
        self,
        *,
        name: str,
        started_at: str,
        wall_s: float,
        peak_rss_mb: float | None,
        exit_code: int,
    ) -> None:
        self.manifest.steps.append(
            StepRecord(
                name=name,
                started_at=started_at,
                wall_s=wall_s,
                peak_rss_mb=peak_rss_mb,
                exit_code=exit_code,
            )
        )

    def warn(self, message: str) -> None:
        if message not in self.manifest.warnings:
            self.manifest.warnings.append(message)

    @contextmanager
    def step(self, name: str) -> Iterator[None]:
        """Time a framework step (peak RSS is the process-lifetime maximum)."""
        started_at = _now()
        t0 = time.perf_counter()
        exit_code = 1
        try:
            yield
            exit_code = 0
        finally:
            self.add_step(
                name=name,
                started_at=started_at,
                wall_s=round(time.perf_counter() - t0, 3),
                peak_rss_mb=_self_peak_rss_mb(),
                exit_code=exit_code,
            )

    def write(self) -> Path:
        self.run_dir.mkdir(parents=True, exist_ok=True)
        data = self.manifest.model_dump(mode="json")
        self.path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
        return self.path
