"""`clinval fidelity --step qc` with samtools and bcftools on the fixture genome."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pysam
import pytest
import yaml
from typer.testing import CliRunner

from insilico_clinval.cli import app
from tests.conftest import FIXTURES

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        shutil.which("samtools") is None or shutil.which("bcftools") is None,
        reason="needs samtools and bcftools",
    ),
]


def _write_bam(path: Path) -> None:
    """Reads tiling chr1:1001-3000 (2 x 50 bp pairs, insert 200) copied from the reference."""
    seq = pysam.FastaFile(str(FIXTURES / "ref.fa")).fetch("chr1")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(seq)}]}
    records = []
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, start in enumerate(range(1000, 2800, 10)):
            mate = start + 150
            for pos, mpos, flag in ((start, mate, 99), (mate, start, 147)):
                a = pysam.AlignedSegment(out.header)
                a.query_name = f"r{i}"
                a.flag = flag
                a.reference_id = 0
                a.reference_start = pos
                a.mapping_quality = 60
                a.cigartuples = [(0, 50)]
                a.query_sequence = seq[pos : pos + 50]
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a.next_reference_id = 0
                a.next_reference_start = mpos
                a.template_length = 200 if flag == 99 else -200
                records.append(a)
        for a in sorted(records, key=lambda r: r.reference_start):
            out.write(a)
    pysam.index(str(path))


def test_qc_step(tmp_path: Path) -> None:
    bam = tmp_path / "sim.bam"
    _write_bam(bam)
    bed = tmp_path / "qc.bed"
    bed.write_text("chr1\t1000\t3000\n")
    data = {
        "name": "qc-test",
        "reference": str(FIXTURES / "ref.fa"),
        "confident_bed": str(FIXTURES / "hg002_mini.bed"),
        "qc_regions": str(bed),
        "datasets": [
            {
                "id": ds,
                "label": ds,
                "alignment": str(bam),
                "vcf": str(FIXTURES / "hg002_mini.vcf.gz"),
                "truth_vcf": str(FIXTURES / "hg002_mini.vcf.gz"),
            }
            for ds in ("simulated", "real")
        ],
        "baseline": "real",
        "happy": {"pass_only": False},
        "threads": 1,
    }
    config = tmp_path / "fidelity.yaml"
    config.write_text(yaml.safe_dump(data))
    out = tmp_path / "out"
    result = CliRunner().invoke(
        app, ["fidelity", "--config", str(config), "--out", str(out), "--step", "qc"]
    )
    assert result.exit_code == 0, result.output
    res = json.loads((out / "fidelity.json").read_text())
    qc = {q["dataset"]: q for q in res["read_qc"]}
    assert qc["simulated"] == {**qc["real"], "dataset": "simulated"}
    sim = qc["simulated"]
    assert sim["reads_total"] == 360
    assert sim["reads_mapped"] == 1.0
    assert sim["error_rate"] == 0.0
    assert sim["insert_size_mean"] == pytest.approx(200, abs=1)
    # 180 pairs x 100 bp over 2000 bp = 9x on average
    assert sim["coverage_mean"] == pytest.approx(9.0, abs=0.1)
    calls = {q["dataset"]: q for q in res["callset_qc"]}
    assert calls["real"]["snvs"] + calls["real"]["indels"] > 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert {t["name"] for t in manifest["tools"]} == {"bcftools", "samtools"}
    assert "alignment:simulated" in {i["role"] for i in manifest["inputs"]}
