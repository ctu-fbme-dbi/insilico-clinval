"""simulate: Mason arguments, coordinate lift, FASTQ merge, planting check."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Any

import pysam
import pytest
from pydantic import ValidationError

from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.enums import Haplotype, PlantingStatus, Zygosity
from insilico_clinval.schemas.sample import SampleVariant
from insilico_clinval.schemas.simulation import SegmentRecord, SimulationProfile
from insilico_clinval.schemas.suite import SimulationConfig
from insilico_clinval.simprep.prepare import TruthRecord
from insilico_clinval.simulate.check import (
    check_variant,
    classify_by_sequence,
    classify_read_allele,
    edit_distance,
    right_shift,
    variant_window,
)
from insilico_clinval.simulate.mason import fragment_count, mason_args, mason_seed
from insilico_clinval.simulate.merge import lift_alignments, merge_fastq_gz

PROFILE: dict[str, Any] = {
    "schema_version": 1,
    "name": "t",
    "technology": "illumina",
    "read_length": 100,
    "fragment": {"model": "normal", "mean": 300, "sd": 30},
    "coverage": 30,
}


def test_profile_validation() -> None:
    SimulationProfile.model_validate(PROFILE)
    with pytest.raises(ValidationError, match="uniform fragment model needs min and max"):
        SimulationProfile.model_validate(
            {**PROFILE, "fragment": {"model": "uniform", "mean": 300, "sd": 0}}
        )
    with pytest.raises(ValidationError, match="greater than or equal to 65536"):
        SimulationProfile.model_validate({**PROFILE, "chunk_size": 500})
    with pytest.raises(ValidationError, match="at least the read length"):
        SimulationProfile.model_validate({**PROFILE, "fragment": {"mean": 80, "sd": 5}})
    with pytest.raises(ValidationError, match="Extra inputs"):
        SimulationProfile.model_validate({**PROFILE, "seed": 3})
    with pytest.raises(ValidationError, match="region mode only"):
        SimulationConfig(profile="p", mode="genome", extra_regions_bed="x.bed")  # type: ignore[arg-type]


def test_fragment_count_and_seed() -> None:
    """N = ceil(coverage/2 * ploidy * L / (2 * read_length))."""
    profile = SimulationProfile.model_validate(PROFILE)
    assert fragment_count(profile, 10_000, 2) == 1500  # 30x over 10 kb with 2 x 100 bp
    assert fragment_count(profile, 10_000, 1) == 750  # haploid: half depth
    assert fragment_count(profile, 1, 1) == 1
    assert mason_seed(7, "s1", "diploid") == mason_seed(7, "s1", "diploid")
    assert mason_seed(7, "s1", "diploid") != mason_seed(7, "s1", "haploid")
    assert 0 <= mason_seed(8, "s1", "diploid") < 2**31 - 1


def test_mason_args() -> None:
    profile = SimulationProfile.model_validate(
        {
            **PROFILE,
            "fragment": {"model": "uniform", "mean": 300, "sd": 0, "min": 250, "max": 350},
            "illumina": {"prob_mismatch": 0.002, "prob_insert": 0.0001},
            "threads": 2,
        }
    )
    args = mason_args(
        profile,
        fasta="g.fa",
        vcf=None,
        fragments=10,
        seed=5,
        read_name_prefix="s.diploid.",
        out_r1="1.fq",
        out_r2="2.fq",
        out_alignment="a.bam",
    )
    text = " ".join(args)
    assert "-iv" not in args
    assert "-ir g.fa -n 10 --seed 5 --num-threads 2 --chunk-size 65536" in text
    assert "--fragment-size-model uniform" in text
    assert "--fragment-min-size 250 --fragment-max-size 350" in text
    assert "--illumina-prob-mismatch 0.002" in text and "--illumina-prob-insert 0.0001" in text
    assert "--illumina-prob-deletion" not in text
    assert args[-8:] == [
        "--read-name-prefix",
        "s.diploid.",
        "-o",
        "1.fq",
        "-or",
        "2.fq",
        "-oa",
        "a.bam",
    ]


def test_merge_fastq_gz_is_deterministic(tmp_path: Path) -> None:
    a, b = tmp_path / "a.fq", tmp_path / "b.fq"
    a.write_text("@r1\nACGT\n+\nIIII\n")
    b.write_text("@r2\nTTTT\n+\nIIII\n")
    out1 = merge_fastq_gz([a, b], tmp_path / "1.fq.gz")
    out2 = merge_fastq_gz([a, b], tmp_path / "sub_dir_name.fq.gz")
    assert out1.read_bytes() == out2.read_bytes()
    assert gzip.decompress(out1.read_bytes()).decode() == a.read_text() + b.read_text()


def test_merge_fastq_gz_does_not_depend_on_threads(tmp_path: Path) -> None:
    # More than one 1 MiB gzip block, split unevenly over two files.
    record = "@r\nACGTTGCAACGTTGCA\n+\nIIIIHHHHGGGGFFFF\n"
    a, b = tmp_path / "a.fq", tmp_path / "b.fq"
    a.write_text(record * 50_000)
    b.write_text(record * 37_123)
    one = merge_fastq_gz([a, b], tmp_path / "t1.fq.gz", threads=1)
    four = merge_fastq_gz([a, b], tmp_path / "t4.fq.gz", threads=4)
    assert one.read_bytes() == four.read_bytes()
    assert gzip.decompress(four.read_bytes()).decode() == a.read_text() + b.read_text()


def _segment_bam(path: Path, segments: list[SegmentRecord]) -> None:
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.4"}, "SQ": [{"SN": s.name, "LN": s.end - s.start + 1} for s in segments]}
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, flag, ref, pos, mate_pos in (
            ("r1", 99, "chr1_101_600", 11, 201),
            ("r1", 147, "chr1_101_600", 201, 11),
        ):
            rec = pysam.AlignedSegment.from_dict(
                {
                    "name": name,
                    "flag": str(flag),
                    "ref_name": ref,
                    "ref_pos": str(pos),
                    "map_quality": "255",
                    "cigar": "10M",
                    "next_ref_name": "=",
                    "next_ref_pos": str(mate_pos),
                    "length": "200",
                    "seq": "ACGTACGTAC",
                    "qual": "IIIIIIIIII",
                    "tags": ["oH:i:2"],
                },
                header,
            )
            out.write(rec)


def test_lift_alignments(tmp_path: Path, reference: Reference) -> None:
    """Contig, position and mate position are lifted back to the genome."""
    segments = [SegmentRecord(name="chr1_101_600", chrom="chr1", start=101, end=600, ploidy=2)]
    src = tmp_path / "seg.bam"
    _segment_bam(src, segments)
    dst = tmp_path / "genome.bam"
    assert lift_alignments([src], segments, reference, dst) == 2
    with pysam.AlignmentFile(str(dst)) as bam:
        assert [s["SN"] for s in bam.header.to_dict()["SQ"]] == ["chr1", "chr2", "chrX", "chrY"]
        recs = list(bam)
    assert [
        (r.reference_name, r.reference_start + 1, r.next_reference_start + 1) for r in recs
    ] == [
        ("chr1", 111, 301),
        ("chr1", 301, 111),
    ]
    assert recs[0].next_reference_name == "chr1" and recs[0].get_tag("oH") == 2
    threaded = tmp_path / "genome_t4.bam"
    lift_alignments([src], segments, reference, threaded, threads=4)
    assert threaded.read_bytes() == dst.read_bytes()


def test_lift_alignments_logs_progress(
    tmp_path: Path, reference: Reference, caplog: pytest.LogCaptureFixture
) -> None:
    segments = [SegmentRecord(name="chr1_101_600", chrom="chr1", start=101, end=600, ploidy=2)]
    src = tmp_path / "seg.bam"
    _segment_bam(src, segments)
    with caplog.at_level(logging.INFO, logger="insilico_clinval"):
        lift_alignments(
            [src], segments, reference, tmp_path / "g.bam", expected=4, progress="s1", log_every=1
        )
        lift_alignments([src], segments, reference, tmp_path / "quiet.bam", log_every=1)
    assert caplog.messages == [
        "s1: lifted 1 / 4 alignments (25%)",
        "s1: lifted 2 / 4 alignments (50%)",
    ]


def test_right_shift_and_classification(reference: Reference) -> None:
    # Fixture homopolymer chr2:9500-9505 = AAAAAA, chr2:9499 = C, chr2:9506 = G.
    assert right_shift(reference, "chr2", 9499, "CA", "C") == 5
    assert right_shift(reference, "chr2", 9499, "C", "CA") == 6
    assert right_shift(reference, "chr2", 9499, "C", "T") == 0
    assert edit_distance("ACGT", "AGT") == 1
    assert classify_read_allele("AAA", "AAA", "AA") == "ref"
    assert classify_read_allele("AAT", "AAA", "TTT") == "ref"
    assert classify_read_allele("GGG", "AAA", "TTT") == "other"


def _planting_bam(
    path: Path, reference: Reference, reads: list[tuple[str, int, str, str, int]]
) -> None:
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": n, "LN": ln} for n, ln in reference.contigs],
        }
    )
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for i, (chrom, start1, cigar, seq, hap) in enumerate(reads):
            rec = pysam.AlignedSegment(header)
            rec.query_name = f"r{i}"
            rec.reference_id = header.get_tid(chrom)
            rec.reference_start = start1 - 1
            rec.cigarstring = cigar
            rec.query_sequence = seq
            rec.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            rec.mapping_quality = 60
            rec.set_tag("oH", hap)
            out.write(rec)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))


def test_check_variant_deletion_in_homopolymer(tmp_path: Path, reference: Reference) -> None:
    """Deletion placed anywhere in the homopolymer counts as ALT; haplotypes separated."""
    start = 9480
    ref_seq = reference.fetch("chr2", start - 1, start - 1 + 40)
    del_at = 9500 - start  # delete the first A of the homopolymer
    alt_seq = ref_seq[:del_at] + ref_seq[del_at + 1 :]
    reads = [
        ("chr2", start, f"{del_at}M1D{39 - del_at}M", alt_seq, 1),  # gap left-aligned
        ("chr2", start, f"{del_at + 5}M1D{34 - del_at}M", alt_seq, 1),  # gap at the right end
        ("chr2", start, "40M", ref_seq, 2),
        ("chr2", start, "40M", ref_seq, 2),
        ("chr2", 9502, "20M", reference.fetch("chr2", 9501, 9521), 2),  # does not span: ignored
    ]
    bam_path = tmp_path / "t.bam"
    _planting_bam(bam_path, reference, reads)
    record = TruthRecord("chr2", 9499, "CA", "C", (1, 0), "v1")
    sv = SampleVariant(variant_id="v1", zygosity=Zygosity.HET, haplotype=Haplotype.ONE)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        result = check_variant(bam, reference, record, sv, min_alt_reads=1)
        assert result.status == PlantingStatus.OK
        assert (result.spanning_reads, result.alt_reads, result.ref_reads) == (4, 2, 2)
        assert [(h.haplotype, h.alt, h.ref) for h in result.by_haplotype] == [(1, 2, 0), (2, 0, 2)]
        assert result.vaf == 0.5 and result.expected_vaf == 0.5
        hom = SampleVariant(variant_id="v1", zygosity=Zygosity.HOM, haplotype=Haplotype.BOTH)
        missing = check_variant(bam, reference, record, hom, min_alt_reads=1)
        assert missing.status == PlantingStatus.MISSING_ON_HAPLOTYPE
        far = TruthRecord("chr1", 5000, reference.fetch("chr1", 4999, 5000), "T", (1, 0), "v2")
        none = check_variant(bam, reference, far, sv, min_alt_reads=1)
        assert none.status == PlantingStatus.NO_READS and none.vaf is None


def test_check_variant_tandem_duplication_scattered_cigar(
    tmp_path: Path, reference: Reference
) -> None:
    """A tandem duplication whose true-alignment CIGAR scatters the inserted bases is still ALT."""
    pos, start = 2000, 1980
    unit = reference.fetch("chr1", pos, pos + 6)  # duplicated right after POS
    ref_seq = reference.fetch("chr1", start - 1, start - 1 + 70)
    k = pos - start + 1
    alt_seq = (ref_seq[:k] + unit + ref_seq[k:])[:70]
    # Mason-like CIGAR: a gap at the left window edge, the inserted bases spread as 1-2 bp gaps
    ref_base = reference.fetch("chr1", pos - 1, pos)
    start0, _ = variant_window(reference, "chr1", pos, ref_base, ref_base + unit)
    lead = start0 - (start - 1)
    ops = [("M", lead), ("D", 1), ("M", 4), ("I", 2), ("M", 3), ("I", 2), ("M", 3), ("I", 3)]
    used = sum(n for op, n in ops if op in "MI")
    cigar = "".join(f"{n}{op}" for op, n in [*ops, ("M", len(alt_seq) - used)])
    reads = [
        ("chr1", start, cigar, alt_seq, 1),
        ("chr1", start, cigar, alt_seq, 1),
        ("chr1", start, "70M", ref_seq, 2),
    ]
    bam_path = tmp_path / "dup.bam"
    _planting_bam(bam_path, reference, reads)
    record = TruthRecord("chr1", pos, ref_base, ref_base + unit, (1, 0), "dup")
    sv = SampleVariant(variant_id="dup", zygosity=Zygosity.HET, haplotype=Haplotype.ONE)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        result = check_variant(bam, reference, record, sv, min_alt_reads=1)
    assert [(h.haplotype, h.alt, h.ref) for h in result.by_haplotype] == [(1, 2, 0), (2, 0, 1)]
    assert result.status == PlantingStatus.OK
    assert classify_by_sequence(alt_seq, "ACGT", "TTTT") is None  # neither window: undecided


def test_check_variant_long_deletion_by_breakpoints(tmp_path: Path, reference: Reference) -> None:
    """A 40 bp deletion: no 100 bp read spans its window, local haplotype sequences decide."""
    pos, k = 3000, 40
    ref_allele = reference.fetch("chr1", pos - 1, pos + k)  # anchor + 40 deleted bases
    alt_allele = ref_allele[0]
    s0 = pos - 1 - 30
    hap_ref = reference.fetch("chr1", s0, s0 + 100)
    after = reference.fetch("chr1", pos + k, pos + k + 70)
    hap_alt = reference.fetch("chr1", s0, pos) + after  # left flank + anchor + bases after
    reads = [
        ("chr1", s0 + 1, "100M", hap_alt[:100], 1),
        ("chr1", s0 + 1, "100M", hap_alt[:100], 1),
        ("chr1", s0 + 1, "100M", hap_ref, 2),
    ]
    bam_path = tmp_path / "long.bam"
    _planting_bam(bam_path, reference, reads)
    record = TruthRecord("chr1", pos, ref_allele, alt_allele, (1, 0), "del40")
    sv = SampleVariant(variant_id="del40", zygosity=Zygosity.HET, haplotype=Haplotype.ONE)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        result = check_variant(bam, reference, record, sv, min_alt_reads=1)
    assert [(h.haplotype, h.alt, h.ref) for h in result.by_haplotype] == [(1, 2, 0), (2, 0, 1)]
    assert result.status == PlantingStatus.OK


def test_check_variant_repeat_deletion_shifted_starts(tmp_path: Path, reference: Reference) -> None:
    """A deletion in a repeat whose reads all start inside the window is decided by k-mers."""
    record = TruthRecord("chr2", 9499, "CA", "C", (1, 0), "v1")  # homopolymer, shift 5
    start0, _ = variant_window(reference, "chr2", 9499, "CA", "C")
    hap_ref = reference.fetch("chr2", 9474, 9534)
    hap_alt = reference.fetch("chr2", 9474, 9499) + reference.fetch("chr2", 9500, 9535)
    # Mason-like: the sequence begins upstream but the alignment start lies inside the window
    reads = [
        ("chr2", start0 + 2, "60M", hap_alt, 1),
        ("chr2", start0 + 2, "60M", hap_alt, 1),
        ("chr2", start0 + 2, "60M", hap_ref, 2),
    ]
    bam_path = tmp_path / "rep.bam"
    _planting_bam(bam_path, reference, reads)
    sv = SampleVariant(variant_id="v1", zygosity=Zygosity.HET, haplotype=Haplotype.ONE)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        result = check_variant(bam, reference, record, sv, min_alt_reads=1)
    assert [(h.haplotype, h.alt, h.ref) for h in result.by_haplotype] == [(1, 2, 0), (2, 0, 1)]
    assert result.status == PlantingStatus.OK
