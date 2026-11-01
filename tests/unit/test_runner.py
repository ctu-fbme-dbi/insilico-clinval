"""Runner abstraction and runtime configuration."""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import pytest
from pydantic import ValidationError
from typer.testing import CliRunner

from insilico_clinval._exec import ToolError
from insilico_clinval._runner import RuntimeConfigError, Toolbox, collapse_mounts, load_runtime
from insilico_clinval.cli import app
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas.runtime import RuntimeConfig
from tests.conftest import EREPO_SUITE, requires_bcftools

DOCKER = {
    "schema_version": 1,
    "default_runner": "local",
    "tools": {"mason_simulator": {"runner": "docker", "image": "mason@sha256:abc"}},
}


def test_runtime_validation_errors() -> None:
    with pytest.raises(ValidationError, match="unknown tool 'masonn'"):
        RuntimeConfig.model_validate({"tools": {"masonn": {"image": "x"}}})
    with pytest.raises(ValidationError, match="needs an image"):
        RuntimeConfig.model_validate({"tools": {"samtools": {"runner": "docker"}}})
    with pytest.raises(ValidationError, match="needs an image"):
        RuntimeConfig.model_validate({"default_runner": "docker", "tools": {"samtools": {}}})
    with pytest.raises(ValidationError, match="Extra inputs"):
        RuntimeConfig.model_validate({"runner": "docker"})


def test_load_runtime(tmp_path: Path) -> None:
    assert load_runtime(None) == RuntimeConfig()
    path = tmp_path / "rt.yaml"
    path.write_text("schema_version: 1\ndefault_runner: dockr\n")
    with pytest.raises(RuntimeConfigError, match="invalid runtime config"):
        load_runtime(path)
    with pytest.raises(RuntimeConfigError, match="not found"):
        load_runtime(tmp_path / "missing.yaml")


def test_collapse_mounts(tmp_path: Path) -> None:
    a = tmp_path / "a"
    nested = a / "b" / "c"
    nested.mkdir(parents=True)
    other = tmp_path / "x"
    other.mkdir()
    assert collapse_mounts([nested, a, other, a]) == sorted([a.resolve(), other.resolve()], key=str)


def test_local_and_docker_commands(tmp_path: Path) -> None:
    box = Toolbox(RuntimeConfig.model_validate(DOCKER))
    assert box.command("samtools", ["index", "x.bam"]) == ["samtools", "index", "x.bam"]
    work = tmp_path / "work"
    (work / "sub").mkdir(parents=True)
    argv = box.command("mason_simulator", ["-n", "5"], mounts=[work / "sub", work], cwd=work)
    resolved = str(work.resolve())
    assert argv[:5] == ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}"]
    assert argv.count("-v") == 1
    assert f"{resolved}:{resolved}" in argv
    assert argv[argv.index("-w") + 1] == resolved
    assert argv[-4:] == ["mason@sha256:abc", "mason_simulator", "-n", "5"]
    with_env = box.command("mason_simulator", ["-n", "5"], cwd=work, env={"TMPDIR": resolved})
    assert with_env[with_env.index("-e") + 1] == f"TMPDIR={resolved}"
    assert with_env.index("-e") < with_env.index("mason@sha256:abc")
    assert box.command("samtools", ["index"], env={"TMPDIR": "/x"}) == ["samtools", "index"]
    assert box.container_digest("mason_simulator") == "sha256:abc"
    assert box.container_digest("samtools") is None
    with pytest.raises(ToolError, match="unknown tool"):
        box.command("gatk", [])


def _fake_tool(directory: Path, name: str, body: str) -> None:
    path = directory / name
    path.write_text(f"#!/bin/sh\n{body}\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _manifest(tmp_path: Path) -> ManifestBuilder:
    return ManifestBuilder(tmp_path / "run", config_path=tmp_path, config_sha256="0" * 64, seed=1)


def test_toolbox_records_version_and_step(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _fake_tool(
        bin_dir,
        "samtools",
        'if [ "$1" = "--version" ]; then echo "samtools 9.9"; echo extra; exit 0; fi\n'
        'if [ "$1" = "fail" ]; then echo "bad input" >&2; exit 3; fi\necho ok',
    )
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    manifest = _manifest(tmp_path)
    box = Toolbox()
    result = box.run("samtools", ["view"], step="samtools view (s1)", manifest=manifest)
    assert result.stdout.strip() == "ok"
    assert [(t.name, t.version, t.container_digest) for t in manifest.manifest.tools] == [
        ("samtools", "samtools 9.9", None)
    ]
    step = manifest.manifest.steps[-1]
    assert step.name == "samtools view (s1)" and step.exit_code == 0
    with pytest.raises(ToolError, match="bad input"):
        box.run("samtools", ["fail"], step="fails", manifest=manifest)
    assert manifest.manifest.steps[-1].exit_code == 3


def test_unpinned_docker_image_is_warned(tmp_path: Path) -> None:
    runtime = RuntimeConfig.model_validate(
        {"tools": {"mason_simulator": {"runner": "docker", "image": "mason:2.4.0"}}}
    )
    box = Toolbox(runtime)
    box._versions["mason_simulator"] = "mason_simulator version: 2.0.9"
    box._digests["mason_simulator"] = "sha256:feed"
    manifest = _manifest(tmp_path)
    box.record("mason_simulator", manifest)
    assert manifest.manifest.tools[0].container_digest == "sha256:feed"
    assert any("not pinned by digest" in w for w in manifest.manifest.warnings)


@pytest.mark.integration
@requires_bcftools
def test_build_truth_with_runtime(tmp_path: Path) -> None:
    """--runtime is accepted and recorded as a manifest input."""
    runtime = tmp_path / "local.yaml"
    runtime.write_text("schema_version: 1\ndefault_runner: local\n")
    out = tmp_path / "run"
    result = CliRunner().invoke(
        app,
        ["build-truth", "--suite", str(EREPO_SUITE), "--out", str(out), "--runtime", str(runtime)],
    )
    assert result.exit_code == 0, result.output
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert "runtime_config" in {i["role"] for i in manifest["inputs"]}
    assert {t["name"] for t in manifest["tools"]} == {"bcftools"}
