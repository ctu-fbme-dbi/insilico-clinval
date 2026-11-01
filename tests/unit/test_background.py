"""Background merge: record filters for assembly-based backgrounds (GIAB v5.0q)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pysam

from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.suite import BackgroundConfig
from insilico_clinval.truth.background import BackgroundMerger, BackgroundRecord

# (chrom, pos, ref, alts, filter, gt, phased); positions on the fixture reference (tests/fixtures).
RECORDS: list[tuple[str, int, str, tuple[str, ...], str | None, tuple[int, ...], bool]] = [
    ("chr1", 100, "T", ("C",), None, (0, 1), True),  # plain het, inside the BED
    ("chr1", 200, "A", ("G", "*"), None, (1, 2), True),  # het SNV + spanning deletion
    ("chr1", 210, "C", ("*",), None, (0, 1), True),  # only a spanning deletion
    ("chr1", 300, "G", ("T",), "GAP1", (0, 1), True),  # filtered
    ("chr1", 400, "A" + "C" * 25, ("A",), None, (1, 1), True),  # 25 bp deletion
    ("chr1", 995, "GACGTACGTA", ("G",), None, (0, 1), True),  # crosses the BED end (1000)
    ("chr1", 1500, "T", ("A",), None, (0, 1), True),  # outside the BED
    ("chrX", 3000, "C", ("T",), None, (1,), False),  # haploid (non-PAR chrX, male)
]


def _write_background(directory: Path) -> tuple[Path, Path]:
    header = pysam.VariantHeader()
    for name, length in (("chr1", 10000), ("chr2", 10000), ("chrX", 8000), ("chrY", 2000)):
        header.contigs.add(name, length=length)
    header.add_line('##FILTER=<ID=GAP1,Description="Uncalled in the first haplotype">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_sample("HG002")
    vcf = directory / "bg.vcf.gz"
    with pysam.VariantFile(str(vcf), "wz", header=header) as out:
        for chrom, pos, ref, alts, filt, gt, phased in RECORDS:
            rec = out.new_record(contig=chrom, start=pos - 1, alleles=(ref, *alts))
            if filt is not None:
                rec.filter.add(filt)
            rec.samples["HG002"]["GT"] = gt
            rec.samples["HG002"].phased = phased
            out.write(rec)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    bed = directory / "bg.bed"
    bed.write_text("chr1\t50\t1000\nchrX\t2700\t3400\n")
    return vcf, bed


def _merge(tmp_path: Path, reference: Reference, **options: Any) -> BackgroundMerger:
    vcf, bed = _write_background(tmp_path)
    config = BackgroundConfig(
        name="GIAB_HG002", vcf=str(vcf), bed=str(bed), release="test", sex="male", **options
    )
    return BackgroundMerger(
        config=config,
        vcf_path=vcf,
        bed_path=bed,
        reference=reference,
        seed=1,
        sample_id="positive",
        planted=[],
        min_distance_bp=0,
    )


def _records(merger: BackgroundMerger) -> list[BackgroundRecord]:
    return [r for chrom in ("chr1", "chrX") for r in merger.records(chrom)]


def test_defaults_keep_the_v421_behaviour(tmp_path: Path, reference: Reference) -> None:
    merger = _merge(tmp_path, reference)
    kept = {(r.chrom, r.pos) for r in _records(merger)}
    assert kept == {
        ("chr1", 100), ("chr1", 300), ("chr1", 400), ("chr1", 995), ("chr1", 1500),
        ("chrX", 3000),
    }  # fmt: skip
    assert merger.skipped == {"non_simple_allele": 2}
    assert merger.spanning_deletion_to_ref == 0


def test_restrict_to_bed(tmp_path: Path, reference: Reference) -> None:
    merger = _merge(tmp_path, reference, restrict_to_bed=True)
    kept = {(r.chrom, r.pos) for r in _records(merger)}
    assert ("chr1", 1500) not in kept and ("chr1", 995) not in kept  # outside / crossing the end
    assert ("chrX", 3000) in kept
    assert merger.skipped["outside_bed"] == 2


def test_require_pass_accepts_pass_and_missing_filter(tmp_path: Path, reference: Reference) -> None:
    merger = _merge(tmp_path, reference, require_pass=True)
    kept = {(r.chrom, r.pos) for r in _records(merger)}
    assert ("chr1", 300) not in kept  # GAP1
    assert ("chr1", 100) in kept  # FILTER missing (as in GIAB v5.0q)
    assert merger.skipped["not_pass"] == 1


def test_max_length(tmp_path: Path, reference: Reference) -> None:
    merger = _merge(tmp_path, reference, max_length=20)
    kept = {(r.chrom, r.pos) for r in _records(merger)}
    assert ("chr1", 400) not in kept
    assert merger.skipped["too_long"] == 1
    merger = _merge(tmp_path, reference, max_length=25)
    assert ("chr1", 400) in {(r.chrom, r.pos) for r in _records(merger)}


def test_spanning_deletion_read_as_ref(tmp_path: Path, reference: Reference) -> None:
    merger = _merge(tmp_path, reference, spanning_deletion="ref")
    by_pos = {(r.chrom, r.pos): r for r in _records(merger)}
    assert by_pos[("chr1", 200)] == BackgroundRecord("chr1", 200, "A", ("G",), (1, 0))
    assert ("chr1", 210) not in by_pos  # only REF and '*' left
    assert merger.skipped == {"no_alt_allele": 1}
    assert merger.spanning_deletion_to_ref == 2


def test_haploid_genotype_kept(tmp_path: Path, reference: Reference) -> None:
    merger = _merge(
        tmp_path, reference, restrict_to_bed=True, require_pass=True, max_length=20,
        spanning_deletion="ref",
    )  # fmt: skip
    chrx = [r for r in _records(merger) if r.chrom == "chrX"]
    assert chrx == [BackgroundRecord("chrX", 3000, "C", ("T",), (1,))]
