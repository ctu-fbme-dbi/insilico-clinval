#!/usr/bin/env python3
"""Download the reference data listed in scripts/data_sources.yaml, with checksums.

    python scripts/fetch_data.py --list
    python scripts/fetch_data.py --only grch38,giab_hg002_v4.2.1
    python scripts/fetch_data.py --only clingen_erepo --allow-unpinned

Every file is downloaded to ``<path>.part``, its SHA-256 is computed while streaming and compared
with the pinned value, and only then renamed to its final path. Existing files are verified, not
downloaded again. A file without a pinned ``sha256`` is fetched only with ``--allow-unpinned``.
Preparation steps (reference decompression + ``.fai``) run after the downloads. Everything done
is written to ``<data-dir>/fetch_report.json``. Exit status: 0 success, 1 failure (download,
checksum mismatch, unpinned file without ``--allow-unpinned``, missing manual dataset).

This is the only component that accesses the network; tests use ``file://`` URLs. Behind a
TLS-inspecting proxy, pass its CA bundle with ``--ca-bundle`` (or set ``SSL_CERT_FILE``).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import ssl
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import pysam
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

DEFAULT_SOURCES = Path(__file__).resolve().parent / "data_sources.yaml"
REPORT_NAME = "fetch_report.json"
_CHUNK = 1 << 20
_TIMEOUT_S = 120


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class FileSpec(_Model):
    url: str
    path: str = Field(description="relative to the data directory; may contain {date}")
    sha256: str | None = None


class PrepareStep(_Model):
    """``decompress_fasta`` (source -> target + .fai), ``run`` (argv; placeholders {data_dir},
    {uid}, {gid}; skipped when ``creates`` exists) or ``remove`` (paths, missing ones ignored)."""

    action: Literal["decompress_fasta", "run", "remove"]
    source: str | None = None
    target: str | None = None
    argv: list[str] = Field(default_factory=list)
    creates: str | None = Field(default=None, description="run: skip if this path exists")
    paths: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _fields(self) -> PrepareStep:
        if self.action == "decompress_fasta" and (self.source is None or self.target is None):
            raise ValueError("decompress_fasta needs source and target")
        if self.action == "run" and not self.argv:
            raise ValueError("run needs argv")
        if self.action == "remove" and not self.paths:
            raise ValueError("remove needs paths")
        return self


class Dataset(_Model):
    description: str
    files: list[FileSpec] = Field(default_factory=list)
    prepare: list[PrepareStep] = Field(default_factory=list)
    manual: str | None = Field(default=None, description="instructions; no download")


class DataSources(_Model):
    schema_version: Literal[1]
    datasets: dict[str, Dataset]


class FetchError(RuntimeError):
    pass


def load_sources(path: Path) -> DataSources:
    return DataSources.model_validate(yaml.safe_load(path.read_text()))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(_CHUNK):
            h.update(chunk)
    return h.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def download(url: str, dest: Path, *, context: ssl.SSLContext | None = None) -> str:
    """Stream ``url`` to ``dest.part``, rename to ``dest``; returns the SHA-256."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    h = hashlib.sha256()
    try:
        with (
            urllib.request.urlopen(url, timeout=_TIMEOUT_S, context=context) as response,
            part.open("wb") as out,
        ):
            while chunk := response.read(_CHUNK):
                h.update(chunk)
                out.write(chunk)
    except OSError as exc:
        part.unlink(missing_ok=True)
        raise FetchError(f"download failed: {url}: {exc}") from exc
    part.rename(dest)
    return h.hexdigest()


def fetch_file(
    spec: FileSpec,
    data_dir: Path,
    *,
    date: str,
    allow_unpinned: bool,
    context: ssl.SSLContext | None = None,
) -> dict[str, object]:
    path = data_dir / spec.path.format(date=date)
    if spec.sha256 is None and not allow_unpinned:
        raise FetchError(
            f"{spec.path}: no pinned sha256 (moving or not yet pinned); rerun with --allow-unpinned"
        )
    action = "verified"
    if path.is_file():
        digest = sha256_file(path)
    else:
        print(f"downloading {spec.url}", file=sys.stderr)
        digest = download(spec.url, path, context=context)
        action = "downloaded"
    if spec.sha256 is not None and digest != spec.sha256:
        bad = path.with_name(path.name + ".sha256-mismatch")
        path.rename(bad)
        raise FetchError(
            f"{path}: SHA-256 mismatch (expected {spec.sha256}, got {digest}); moved to {bad.name}"
        )
    if spec.sha256 is None:
        print(f"unpinned: {path.relative_to(data_dir)} sha256={digest}", file=sys.stderr)
    return {
        "url": spec.url,
        "path": str(path.relative_to(data_dir)),
        "size": path.stat().st_size,
        "sha256": digest,
        "pinned": spec.sha256 is not None,
        "action": action,
        "time": _now(),
    }


def decompress_fasta(source: Path, target: Path) -> None:
    """gunzip a FASTA (plain gzip or bgzip) and index it with pysam (samtools faidx)."""
    if not (target.is_file() and Path(f"{target}.fai").is_file()):
        part = target.with_name(target.name + ".part")
        with gzip.open(source, "rb") as src, part.open("wb") as dst:
            shutil.copyfileobj(src, dst, _CHUNK)
        part.rename(target)
        pysam.faidx(str(target))


def run_prepare(step: PrepareStep, data_dir: Path) -> dict[str, object]:
    if step.action == "decompress_fasta":
        assert step.source is not None and step.target is not None
        target = data_dir / step.target
        decompress_fasta(data_dir / step.source, target)
        return {
            "action": step.action,
            "source": step.source,
            "target": step.target,
            "sha256": sha256_file(target),
            "time": _now(),
        }
    if step.action == "remove":
        removed = []
        for rel in step.paths:
            path = data_dir / rel
            if path.is_file():
                path.unlink()
                removed.append(rel)
        return {"action": "remove", "removed": removed, "time": _now()}
    values = {"data_dir": str(data_dir.resolve()), "uid": str(os.getuid()), "gid": str(os.getgid())}
    argv = [arg.format_map(values) for arg in step.argv]
    if step.creates is not None and (data_dir / step.creates).exists():
        return {"action": "run", "argv": argv, "skipped": f"{step.creates} exists", "time": _now()}
    print("running " + " ".join(argv), file=sys.stderr)
    result = subprocess.run(argv, check=False)
    if result.returncode != 0:
        raise FetchError(f"{argv[0]} exited with status {result.returncode}")
    if step.creates is not None and not (data_dir / step.creates).exists():
        # Some tools (e.g. Nirvana's Downloader on a usage error) exit 0 without output.
        raise FetchError(f"{argv[0]} exited 0 but did not create {step.creates}")
    return {"action": "run", "argv": argv, "time": _now()}


def fetch(
    sources: DataSources,
    data_dir: Path,
    *,
    only: list[str] | None,
    allow_unpinned: bool,
    date: str,
    context: ssl.SSLContext | None = None,
) -> tuple[dict[str, object], list[str]]:
    names = only if only else list(sources.datasets)
    unknown = [n for n in names if n not in sources.datasets]
    if unknown:
        raise FetchError(f"unknown datasets: {', '.join(unknown)}")
    report: dict[str, object] = {}
    errors: list[str] = []
    for name in names:
        dataset = sources.datasets[name]
        if dataset.manual is not None:
            errors.append(f"{name}: manual dataset, not downloaded: {dataset.manual.strip()}")
            report[name] = {"manual": dataset.manual.strip()}
            continue
        entry: dict[str, object] = {"files": [], "prepare": []}
        try:
            entry["files"] = [
                fetch_file(f, data_dir, date=date, allow_unpinned=allow_unpinned, context=context)
                for f in dataset.files
            ]
            entry["prepare"] = [run_prepare(step, data_dir) for step in dataset.prepare]
        except FetchError as exc:
            errors.append(f"{name}: {exc}")
            entry["error"] = str(exc)
        report[name] = entry
    return report, errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--only", help="comma-separated dataset names (default: all)")
    parser.add_argument("--allow-unpinned", action="store_true")
    parser.add_argument("--date", help="value of {date} in paths (default: today, UTC)")
    parser.add_argument("--ca-bundle", type=Path, help="CA certificates (PEM) for HTTPS")
    parser.add_argument("--list", action="store_true", help="list datasets and exit")
    args = parser.parse_args(argv)
    sources = load_sources(args.sources)
    if args.list:
        for name, ds in sources.datasets.items():
            kind = "manual" if ds.manual else f"{len(ds.files)} files"
            pinned = all(f.sha256 for f in ds.files)
            state = "-" if ds.manual else ("pinned" if pinned else "unpinned")
            print(f"{name}\t{kind}\t{state}\t{ds.description.strip()}")
        return 0
    only = [s.strip() for s in args.only.split(",")] if args.only else None
    date = args.date or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    data_dir: Path = args.data_dir
    data_dir.mkdir(parents=True, exist_ok=True)
    context = ssl.create_default_context(cafile=str(args.ca_bundle)) if args.ca_bundle else None
    try:
        report, errors = fetch(
            sources,
            data_dir,
            only=only,
            allow_unpinned=args.allow_unpinned,
            date=date,
            context=context,
        )
    except FetchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    report_path = data_dir / REPORT_NAME
    previous = json.loads(report_path.read_text()) if report_path.is_file() else {}
    previous.update(report)
    report_path.write_text(json.dumps(previous, indent=2, sort_keys=True) + "\n")
    for e in errors:
        print(f"error: {e}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
