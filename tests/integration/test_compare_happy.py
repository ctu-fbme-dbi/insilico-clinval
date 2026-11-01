"""hap.py in its image on the fixture genome. Opt-in:

CLINVAL_TEST_HAPPY_IMAGE=clinval-happy:0.3.15 pytest tests/integration/test_compare_happy.py
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pysam
import pytest

from insilico_clinval._runner import Toolbox
from insilico_clinval.compare.evidence import happy_calls, normalized, pick
from insilico_clinval.compare.run import HappyJob, run_happy
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.runtime import RuntimeConfig, ToolRuntime
from tests.conftest import FIXTURES

IMAGE = os.environ.get("CLINVAL_TEST_HAPPY_IMAGE")
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        IMAGE is None or shutil.which("docker") is None, reason="needs hap.py image"
    ),
]


def _vcf(path: Path, rows: list[str]) -> Path:
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=10000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS",
        *rows,
    ]
    plain = path.with_suffix("")
    plain.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return path


def test_happy_decisions(tmp_path: Path) -> None:
    reference = Reference(FIXTURES / "ref.fa")
    shutil.copy(FIXTURES / "ref.fa", tmp_path / "ref.fa")
    shutil.copy(FIXTURES / "ref.fa.fai", tmp_path / "ref.fa.fai")
    b1, b2 = reference.fetch("chr1", 1000, 1001), reference.fetch("chr1", 3000, 3001)
    a1, a2 = ("A" if b1 != "A" else "C"), ("A" if b2 != "A" else "C")
    truth = _vcf(tmp_path / "truth.vcf.gz", [
        f"chr1\t1001\t.\t{b1}\t{a1}\t.\tPASS\t.\tGT\t0|1",
        f"chr1\t3001\t.\t{b2}\t{a2}\t.\tPASS\t.\tGT\t0|1",
    ])  # fmt: skip
    query = _vcf(tmp_path / "query.vcf.gz", [
        f"chr1\t1001\t.\t{b1}\t{a1}\t50\tPASS\t.\tGT\t0/1",
        f"chr1\t3001\t.\t{b2}\t{a2}\t50\tPASS\t.\tGT\t1/1",
    ])  # fmt: skip
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t0\t10000\n")
    strata = tmp_path / "strata.tsv"
    strata.write_text(f"planted\t{regions}\n")
    assert IMAGE is not None
    toolbox = Toolbox(
        RuntimeConfig(default_runner="docker", tools={"hap.py": ToolRuntime(image=IMAGE)})
    )
    job = HappyJob(truth, query, tmp_path / "ref.fa", regions, strata, tmp_path / "h" / "out", 1)
    run_happy(job, toolbox, None, "hap.py")
    with pysam.VariantFile(str(tmp_path / "h" / "out.vcf.gz")) as vcf:
        ok = pick(happy_calls(vcf, reference, "TRUTH", "chr1", 990, 1010),
                  normalized(reference, ("chr1", 1001, b1, a1)))[0]  # fmt: skip
        gt = pick(happy_calls(vcf, reference, "TRUTH", "chr1", 2990, 3010),
                  normalized(reference, ("chr1", 3001, b2, a2)))[0]  # fmt: skip
    assert ok is not None and ok.decision == "TP"
    assert gt is not None and (gt.decision, gt.kind) == ("FN", "am")
