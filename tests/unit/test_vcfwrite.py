"""Truth VCF content before normalisation."""

from __future__ import annotations

import logging
from pathlib import Path

import pysam
import pytest

from insilico_clinval.config import LoadedSuite
from insilico_clinval.genome.reference import Reference
from insilico_clinval.truth.genotype import assign_genotypes
from insilico_clinval.truth.vcfwrite import write_truth_vcf
from insilico_clinval.variants.select import select_variants, variants_per_sample


def _write(loaded_suite: LoadedSuite, reference: Reference, tmp_path: Path, sample: int) -> Path:
    result = select_variants(loaded_suite, reference)
    cfg = loaded_suite.config
    s = cfg.samples[sample]
    genotypes, _ = assign_genotypes(cfg, s, variants_per_sample(result)[s.sample_id])
    path, _ = write_truth_vcf(tmp_path / f"{s.sample_id}.vcf", s.sample_id, genotypes, reference)
    return path


def test_header_contigs_once_and_format_gt_ps(
    loaded_suite: LoadedSuite, reference: Reference, tmp_path: Path
) -> None:
    path = _write(loaded_suite, reference, tmp_path, 0)
    text = path.read_text()
    contig_lines = [line for line in text.splitlines() if line.startswith("##contig")]
    assert len(contig_lines) == len(reference.contigs) == 4
    assert len({line for line in contig_lines}) == 4
    fmt_lines = [line for line in text.splitlines() if line.startswith("##FORMAT")]
    assert [line.split(",")[0] for line in fmt_lines] == ["##FORMAT=<ID=GT", "##FORMAT=<ID=PS"]
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            assert set(rec.format.keys()) <= {"GT", "PS"}
            assert rec.info["TV_ID"].endswith(f"{rec.pos}:{rec.ref}:{rec.alts[0]}")
            assert "TV_CATEGORY" in rec.info and "TV_RELEASE" in rec.info


def test_genotype_strings(loaded_suite: LoadedSuite, reference: Reference, tmp_path: Path) -> None:
    pos_path = _write(loaded_suite, reference, tmp_path, 0)
    neg_path = _write(loaded_suite, reference, tmp_path, 1)
    body = [line.split("\t") for line in pos_path.read_text().splitlines() if line[0] != "#"]
    gts = {(f[0], int(f[1])): (f[8], f[9]) for f in body}
    assert gts[("chr2", 5100)] == ("GT", "1|1")  # synthetic hom
    fmt, value = gts[("chr1", 200)]  # compound het with the second GENEA variant
    assert fmt == "GT:PS" and value.endswith(":200") and value[:3] in ("0|1", "1|0")
    het = [v for f, v in gts.values() if f == "GT" and v != "1|1"]
    assert het and all(v in ("0|1", "1|0") for v in het)
    neg = [line.split("\t") for line in neg_path.read_text().splitlines() if line[0] != "#"]
    chrx = [f[9] for f in neg if f[0] == "chrX"]
    assert chrx == ["1", "1"]  # hemizygous in the male negative sample


def test_write_logs_progress(
    loaded_suite: LoadedSuite,
    reference: Reference,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    result = select_variants(loaded_suite, reference)
    cfg = loaded_suite.config
    s = cfg.samples[0]
    genotypes, _ = assign_genotypes(cfg, s, variants_per_sample(result)[s.sample_id])
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="insilico_clinval"):
        _, n = write_truth_vcf(tmp_path / "p.vcf", s.sample_id, genotypes, reference, log_every=1)
    assert n == len(genotypes) > 0
    assert len(caplog.messages) == n
    assert caplog.messages[0].startswith(f"{s.sample_id}: truth VCF: 1 records written (at chr")
