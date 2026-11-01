"""scripts/fetch_data.py with file:// URLs only."""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "fetch_data.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("fetch_data", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["fetch_data"] = module
    spec.loader.exec_module(module)
    return module


fetch_data = _load()


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def remote(tmp_path: Path) -> dict[str, Any]:
    """A fake 'remote' directory with a BED file and a gzipped FASTA."""
    src = tmp_path / "remote"
    src.mkdir()
    bed = b"chr1\t0\t100\n"
    (src / "a.bed").write_bytes(bed)
    fasta = b">chr1\nACGTACGTAC\n>chr2\nGGGG\n"
    (src / "ref.fa.gz").write_bytes(gzip.compress(fasta))
    return {"dir": src, "bed_sha": _sha(bed), "fasta": fasta}


def _sources(tmp_path: Path, datasets: dict[str, Any]) -> Path:
    path = tmp_path / "sources.yaml"
    path.write_text(yaml.safe_dump({"schema_version": 1, "datasets": datasets}))
    return path


def _bed_dataset(remote: dict[str, Any], sha: str | None) -> dict[str, Any]:
    url = (remote["dir"] / "a.bed").as_uri()
    return {"description": "bed", "files": [{"url": url, "path": "beds/a.bed", "sha256": sha}]}


def test_verified_download_and_idempotent_rerun(tmp_path: Path, remote: dict[str, Any]) -> None:
    sources = _sources(tmp_path, {"bed": _bed_dataset(remote, remote["bed_sha"])})
    data = tmp_path / "data"
    args = ["--sources", str(sources), "--data-dir", str(data)]
    assert fetch_data.main(args) == 0
    assert (data / "beds/a.bed").read_bytes() == b"chr1\t0\t100\n"
    report = json.loads((data / "fetch_report.json").read_text())
    entry = report["bed"]["files"][0]
    assert (
        entry["action"] == "downloaded" and entry["pinned"] and entry["sha256"] == remote["bed_sha"]
    )
    assert fetch_data.main(args) == 0
    report = json.loads((data / "fetch_report.json").read_text())
    assert report["bed"]["files"][0]["action"] == "verified"


def test_checksum_mismatch_fails_without_target(tmp_path: Path, remote: dict[str, Any]) -> None:
    sources = _sources(tmp_path, {"bed": _bed_dataset(remote, "0" * 64)})
    data = tmp_path / "data"
    assert fetch_data.main(["--sources", str(sources), "--data-dir", str(data)]) == 1
    assert not (data / "beds/a.bed").exists()
    assert (data / "beds/a.bed.sha256-mismatch").exists()
    assert not list(data.glob("beds/*.part"))


def test_unpinned_needs_flag(tmp_path: Path, remote: dict[str, Any]) -> None:
    sources = _sources(tmp_path, {"bed": _bed_dataset(remote, None)})
    data = tmp_path / "data"
    args = ["--sources", str(sources), "--data-dir", str(data)]
    assert fetch_data.main(args) == 1
    assert not (data / "beds/a.bed").exists()
    assert fetch_data.main([*args, "--allow-unpinned"]) == 0
    report = json.loads((data / "fetch_report.json").read_text())
    assert report["bed"]["files"][0]["pinned"] is False


def test_date_placeholder_and_failed_download(tmp_path: Path, remote: dict[str, Any]) -> None:
    url = (remote["dir"] / "a.bed").as_uri()
    missing = (remote["dir"] / "missing.tsv").as_uri()
    sources = _sources(
        tmp_path,
        {
            "dated": {"description": "d", "files": [{"url": url, "path": "x/e_{date}.tsv"}]},
            "gone": {"description": "g", "files": [{"url": missing, "path": "x/gone.tsv"}]},
        },
    )
    data = tmp_path / "data"
    base = ["--sources", str(sources), "--data-dir", str(data), "--allow-unpinned"]
    assert fetch_data.main([*base, "--only", "dated", "--date", "2026-10-07"]) == 0
    assert (data / "x/e_2026-10-07.tsv").is_file()
    assert fetch_data.main([*base, "--only", "gone"]) == 1
    assert not list((data / "x").glob("gone*"))
    assert fetch_data.main([*base, "--only", "nope"]) == 1


def test_prepare_reference(tmp_path: Path, remote: dict[str, Any]) -> None:
    data_bytes = (remote["dir"] / "ref.fa.gz").read_bytes()
    dataset = {
        "description": "ref",
        "files": [
            {
                "url": (remote["dir"] / "ref.fa.gz").as_uri(),
                "path": "reference/ref.fa.gz",
                "sha256": _sha(data_bytes),
            }
        ],
        "prepare": [
            {
                "action": "decompress_fasta",
                "source": "reference/ref.fa.gz",
                "target": "reference/ref.fa",
            }
        ],
    }
    data = tmp_path / "data"
    sources = _sources(tmp_path, {"ref": dataset})
    assert fetch_data.main(["--sources", str(sources), "--data-dir", str(data)]) == 0
    assert (data / "reference/ref.fa").read_bytes() == remote["fasta"]
    fai = (data / "reference/ref.fa.fai").read_text().splitlines()
    assert [line.split("\t")[:2] for line in fai] == [["chr1", "10"], ["chr2", "4"]]


def test_manual_dataset_is_reported(tmp_path: Path) -> None:
    sources = _sources(tmp_path, {"old": {"description": "o", "manual": "ask the authors"}})
    data = tmp_path / "data"
    assert fetch_data.main(["--sources", str(sources), "--data-dir", str(data)]) == 1
    assert json.loads((data / "fetch_report.json").read_text())["old"] == {
        "manual": "ask the authors"
    }


def test_repository_sources_file_is_valid() -> None:
    sources = fetch_data.load_sources(fetch_data.DEFAULT_SOURCES)
    assert {"grch38", "giab_hg002_v4.2.1", "clingen_erepo"} <= set(sources.datasets)
    for ds in sources.datasets.values():
        for f in ds.files:
            assert f.url.startswith("https://")


def test_run_step_must_create_its_output(tmp_path: Path) -> None:
    quiet = [sys.executable, "-c", "pass"]  # exits 0 and writes nothing
    sources = _sources(
        tmp_path,
        {
            "tool": {
                "description": "t",
                "prepare": [{"action": "run", "argv": quiet, "creates": "x"}],
            }
        },
    )
    data = tmp_path / "data"
    assert fetch_data.main(["--sources", str(sources), "--data-dir", str(data)]) == 1
    assert (
        "did not create x" in json.loads((data / "fetch_report.json").read_text())["tool"]["error"]
    )
