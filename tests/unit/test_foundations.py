"""rng, _exec, provenance, alleles, BED intervals."""

from __future__ import annotations

import gzip
import json
import sys
from importlib.metadata import version
from pathlib import Path

import pysam
import pytest

from insilico_clinval._exec import ToolError, run
from insilico_clinval.genome.alleles import classify, event_length, indel_length
from insilico_clinval.genome.bed import IntervalSet
from insilico_clinval.provenance import ManifestBuilder, ManifestError, sha256_file
from insilico_clinval.rng import derive_choice, derive_int
from insilico_clinval.schemas.enums import VariantType


def test_derive_int_is_pure_and_seed_dependent() -> None:
    assert derive_int(1, "rank", "x") == derive_int(1, "rank", "x")
    assert derive_int(1, "rank", "x") != derive_int(2, "rank", "x")
    assert derive_int(1, "rank", "x") != derive_int(1, "haplotype", "x")
    assert 0 <= derive_choice(1, "p", "k", 2) < 2


def test_run_captures_output_and_resources() -> None:
    result = run([sys.executable, "-c", "import sys; print('out'); print('err', file=sys.stderr)"])
    assert result.returncode == 0
    assert result.stdout.strip() == "out"
    assert result.stderr.strip() == "err"
    assert result.wall_s >= 0
    assert result.peak_rss_mb is not None and result.peak_rss_mb > 0


def test_run_raises_on_failure_with_stderr() -> None:
    with pytest.raises(ToolError, match="boom"):
        run([sys.executable, "-c", "import sys; sys.stderr.write('boom'); sys.exit(3)"])
    result = run([sys.executable, "-c", "import sys; sys.exit(3)"], check=False)
    assert result.returncode == 3


def test_run_writes_stderr_file_also_on_failure(tmp_path: Path) -> None:
    log_path = tmp_path / "tool.log"
    code = "import sys; sys.stderr.write('progress ...'); sys.exit(139)"
    result = run([sys.executable, "-c", code], check=False, stderr_path=log_path)
    assert result.returncode == 139
    assert result.stderr == "progress ..." == log_path.read_text()
    with pytest.raises(ToolError, match="boom"):
        run(
            [sys.executable, "-c", "import sys; sys.stderr.write('boom'); sys.exit(1)"],
            stderr_path=log_path,
        )
    assert log_path.read_text() == "boom"


def test_run_passes_extra_env() -> None:
    code = "import os; print(os.environ['TMPDIR'], 'PATH' in os.environ)"
    result = run([sys.executable, "-c", code], env={"TMPDIR": "/scratch/x"})
    assert result.stdout.strip() == "/scratch/x True"


def test_run_missing_tool() -> None:
    with pytest.raises(ToolError, match="not found"):
        run(["definitely-not-a-tool-xyz"])


def test_manifest_extend_and_config_guard(tmp_path: Path) -> None:
    f = tmp_path / "input.txt"
    f.write_text("hello\n")
    m = ManifestBuilder(tmp_path, config_path=f, config_sha256="a" * 64, seed=1)
    m.add_input("x", f)
    m.add_tool("bcftools", "bcftools 1.0")
    with m.step("demo"):
        pass
    m.write()
    m2 = ManifestBuilder(tmp_path, config_path=f, config_sha256="a" * 64, seed=1)
    assert [s.name for s in m2.manifest.steps] == ["demo"]
    assert m2.manifest.inputs[0].sha256 == sha256_file(f)
    with pytest.raises(ManifestError):
        ManifestBuilder(tmp_path, config_path=f, config_sha256="b" * 64, seed=1)


def test_manifest_records_python_packages(tmp_path: Path) -> None:
    f = tmp_path / "input.txt"
    f.write_text("hello\n")
    m = ManifestBuilder(tmp_path, config_path=f, config_sha256="a" * 64, seed=1)
    assert m.manifest.framework.python_packages == {
        "isal": version("isal"),
        "pysam": version("pysam"),
    }
    m.write()
    # A manifest written before python_packages existed, or by another isal version.
    data = json.loads(m.path.read_text())
    data["framework"]["python_packages"] = {"isal": "0.0.1"}
    m.path.write_text(json.dumps(data))
    m2 = ManifestBuilder(tmp_path, config_path=f, config_sha256="a" * 64, seed=1)
    assert m2.manifest.framework.python_packages == {"isal": "0.0.1", "pysam": version("pysam")}
    assert m2.manifest.warnings == [
        f"isal {version('isal')} differs from 0.0.1 recorded earlier in this run"
    ]
    del data["framework"]["python_packages"]
    m.path.write_text(json.dumps(data))
    m3 = ManifestBuilder(tmp_path, config_path=f, config_sha256="a" * 64, seed=1)
    assert set(m3.manifest.framework.python_packages) == {"isal", "pysam"}
    assert m3.manifest.warnings == []


@pytest.mark.parametrize(
    ("ref", "alt", "vtype", "ilen", "elen"),
    [
        ("A", "G", VariantType.SNV, 0, 1),
        ("AC", "GT", VariantType.MNV, 0, 2),
        ("A", "ACGT", VariantType.INSERTION, 3, 3),
        ("ACGT", "A", VariantType.DELETION, -3, 3),
        ("GAG", "CA", VariantType.DELINS, -1, 3),
        ("GAG", "GT", VariantType.DELINS, -1, 2),
    ],
)
def test_classify(ref: str, alt: str, vtype: VariantType, ilen: int, elen: int) -> None:
    assert classify(ref, alt) == vtype
    assert indel_length(ref, alt) == ilen
    assert event_length(ref, alt) == elen


def _bed_rows() -> str:
    return "chr1\t10\t20\nchr1\t15\t30\nchr1\t100\t110\n"


def test_interval_set_plain_and_tabix(tmp_path: Path) -> None:
    plain = tmp_path / "a.bed"
    plain.write_text(_bed_rows())
    gz_noindex = tmp_path / "b.bed.gz"
    with gzip.open(gz_noindex, "wt") as fh:
        fh.write(_bed_rows())
    tabixed = tmp_path / "c.bed"
    tabixed.write_text(_bed_rows())
    tabix_gz = Path(pysam.tabix_index(str(tabixed), preset="bed", force=True))
    for path in (plain, gz_noindex, tabix_gz):
        bed = IntervalSet(path)
        assert bed.overlaps("chr1", 29, 31)  # merged [10, 30)
        assert bed.contains("chr1", 12, 28)
        assert not bed.contains("chr1", 25, 35)
        assert not bed.overlaps("chr1", 30, 100)
        assert not bed.overlaps("chr2", 0, 1000)
