"""simprep: decomposition, conflicts, segments, Mason VCF."""

from __future__ import annotations

import random
from collections.abc import Iterator
from pathlib import Path

import pysam
import pytest

from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.enums import Sex, SimulationMode
from insilico_clinval.schemas.simulation import SegmentRecord
from insilico_clinval.schemas.suite import Region, SimulationConfig
from insilico_clinval.simprep import SimPrepError
from insilico_clinval.simprep.events import Part, apply_parts, conflicts, decompose
from insilico_clinval.simprep.mason_vcf import MasonVcfWriter, validate_mason_vcf
from insilico_clinval.simprep.prepare import batch_segments, prepare_sample, read_truth
from insilico_clinval.simprep.segments import (
    SegmentIndex,
    build_segments,
    genome_intervals,
    region_intervals,
)

PARS = (Region(chrom="chrX", start=1, end=1000), Region(chrom="chrY", start=1, end=500))


# -- decomposition ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("ref", "alt", "expected"),
    [
        ("A", "G", [Part(10, "A", "G")]),
        ("A", "ATTG", [Part(10, "A", "ATTG")]),
        ("ACG", "A", [Part(10, "ACG", "A")]),
        ("AC", "GT", [Part(10, "A", "G"), Part(11, "C", "T")]),
        ("GAG", "TT", [Part(10, "G", "T"), Part(11, "A", "T"), Part(11, "AG", "A")]),
        ("A", "GT", [Part(10, "A", "AT"), Part(10, "A", "G")]),
        ("ACT", "AGTCC", [Part(11, "C", "G"), Part(12, "T", "TCC")]),
    ],
)
def test_decompose_cases(ref: str, alt: str, expected: list[Part]) -> None:
    assert decompose(10, ref, alt) == expected


def test_decompose_reproduces_alt_haplotype() -> None:
    """Random MNV/delins/indels; applied parts equal the ALT haplotype; no complex parts."""
    rng = random.Random(3)
    for _ in range(2000):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        pos = rng.randint(5, 30)
        ref = seq[pos - 1 : pos - 1 + rng.randint(1, 8)]
        alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 8)))
        if ref == alt:
            continue
        parts = decompose(pos, ref, alt)
        assert all(len(p.ref) == 1 or len(p.alt) == 1 for p in parts)
        assert all(p.is_indel is False or p.ref[0] == p.alt[0] for p in parts)
        expected = seq[: pos - 1] + alt + seq[pos - 1 + len(ref) :]
        assert apply_parts(seq, parts) == expected, (pos, ref, alt, parts)


def test_conflicts() -> None:
    snv = Part(10, "A", "G")  # span [9, 10)
    assert conflicts(snv, Part(10, "A", "C"))
    assert not conflicts(snv, Part(11, "C", "T"))
    deletion = Part(20, "ACG", "A")  # span [20, 22)
    assert conflicts(deletion, Part(21, "C", "T"))
    assert not conflicts(deletion, Part(20, "A", "T"))  # SNV at the anchor
    assert not conflicts(deletion, Part(23, "T", "A"))
    insertion = Part(30, "A", "AT")  # point 30
    assert not conflicts(insertion, Part(30, "A", "G"))  # SNV at the anchor
    assert not conflicts(insertion, Part(31, "C", "G"))  # SNV right after the insertion point
    assert conflicts(insertion, Part(30, "A", "AG"))  # two insertions at one point
    assert conflicts(insertion, Part(30, "AC", "A"))  # deletion starting at the point
    assert conflicts(Part(29, "AAC", "A"), insertion)  # insertion inside a deletion


# -- segments -------------------------------------------------------------------------------


def test_region_segments_clip_merge_and_ploidy(reference: Reference) -> None:
    """Flanks clipped to contigs, overlapping regions merged, PAR split in males."""
    variants = [
        ("chr1", 100, "A"),
        ("chr1", 900, "ACG"),
        ("chr2", 9950, "T"),
        ("chrX", 950, "G"),
        ("chrY", 1500, "C"),
    ]
    intervals = region_intervals(variants, reference, 600)
    male = build_segments(intervals, reference, Sex.MALE, PARS)
    assert [(s.chrom, s.start, s.end, s.ploidy) for s in male] == [
        ("chr1", 1, 1502, 2),
        ("chr2", 9350, 10000, 2),
        ("chrX", 350, 1000, 2),
        ("chrX", 1001, 1550, 1),
        ("chrY", 900, 2000, 1),
    ]
    assert male[0].name == "chr1_1_1502"
    female = build_segments(intervals, reference, Sex.FEMALE, PARS)
    assert [s.chrom for s in female] == ["chr1", "chr2", "chrX"]
    assert female[-1].ploidy == 2


def test_genome_segments(reference: Reference) -> None:
    intervals = genome_intervals(reference, None)
    male = build_segments(intervals, reference, Sex.MALE, PARS)
    assert [(s.chrom, s.start, s.end, s.ploidy) for s in male] == [
        ("chr1", 1, 10000, 2),
        ("chr2", 1, 10000, 2),
        ("chrX", 1, 1000, 2),
        ("chrX", 1001, 8000, 1),
        ("chrY", 501, 2000, 1),  # chrY PAR left out (represented by the chrX PAR)
    ]
    with pytest.raises(ValueError, match="not in the reference"):
        genome_intervals(reference, ["chr3"])


def test_segment_index() -> None:
    segs = [
        SegmentRecord(name="a", chrom="chr1", start=10, end=20, ploidy=2),
        SegmentRecord(name="b", chrom="chr1", start=40, end=50, ploidy=2),
    ]
    index = SegmentIndex(segs)
    assert index.overlapping("chr1", 15, 15) == segs[0]
    assert index.overlapping("chr1", 30, 35) is None
    assert index.overlapping("chr1", 18, 42) in segs
    assert index.contains("chr1", 40, 50) and not index.contains("chr1", 18, 22)
    assert index.overlapping("chr2", 15, 15) is None


# -- truth VCF -> Mason inputs --------------------------------------------------------------

Rec = tuple[str, int, str, str, tuple[int, ...], str | None]


def vid_of(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"curated:{chrom}:{pos}:{ref}:{alt}"


def planted_ids(records: list[Rec]) -> list[str]:
    return [vid_of(c, p, r, a) for c, p, r, a, _, v in records if v is not None]


def _truth_vcf(path: Path, reference: Reference, records: list[Rec]) -> Path:
    header = pysam.VariantHeader()
    for name, length in reference.contigs:
        header.contigs.add(name, length=length)
    header.add_line('##INFO=<ID=TV_ID,Number=1,Type=String,Description="id">')
    header.add_line('##INFO=<ID=TV_BG,Number=0,Type=Flag,Description="bg">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    header.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">')
    header.add_sample("s1")
    plain = path.with_suffix("")
    with pysam.VariantFile(str(plain), "w", header=header) as out:
        for chrom, pos, ref, alt, gt, vid in records:
            r = out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt))
            if vid is None:
                r.info["TV_BG"] = True
            else:
                r.info["TV_ID"] = vid_of(chrom, pos, ref, alt)
            r.samples["s1"]["GT"] = gt
            r.samples["s1"].phased = len(gt) == 2
            if vid is not None and len(gt) == 2:
                r.samples["s1"]["PS"] = pos
            out.write(r)
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def _ref(reference: Reference, chrom: str, pos: int, n: int = 1) -> str:
    return reference.fetch(chrom, pos - 1, pos - 1 + n)


def _other(base: str) -> str:
    return "A" if base != "A" else "C"


@pytest.fixture
def records(reference: Reference) -> Iterator[list[Rec]]:
    r = reference
    yield [
        (
            "chr1",
            1000,
            _ref(r, "chr1", 1000, 3),
            _other(_ref(r, "chr1", 1000)) + "T",
            (0, 1),
            "p:del",
        ),
        ("chr1", 1002, _ref(r, "chr1", 1002), _other(_ref(r, "chr1", 1002)), (0, 1), None),
        ("chr1", 1001, _ref(r, "chr1", 1001), _other(_ref(r, "chr1", 1001)), (1, 0), None),
        ("chr1", 1300, _ref(r, "chr1", 1300), _ref(r, "chr1", 1300) + "GG", (1, 1), None),
        ("chr1", 3000, _ref(r, "chr1", 3000), _other(_ref(r, "chr1", 3000)), (1, 1), None),
        ("chr2", 4000, _ref(r, "chr2", 4000), _other(_ref(r, "chr2", 4000)), (1, 0), "p:snv"),
        ("chrX", 3000, _ref(r, "chrX", 3000), _other(_ref(r, "chrX", 3000)), (1,), "p:hemi"),
        ("chrX", 3100, _ref(r, "chrX", 3100), _other(_ref(r, "chrX", 3100)), (0, 1), None),
        ("chrX", 3200, _ref(r, "chrX", 3200), _other(_ref(r, "chrX", 3200)), (1, 1), None),
    ]


def _sorted(records: list[Rec], reference: Reference) -> list[Rec]:
    return sorted(records, key=lambda x: (reference.order[x[0]], x[1]))


def _prepare(tmp_path: Path, reference: Reference, records: list[Rec], **kw: object):  # type: ignore[no-untyped-def]
    truth = _truth_vcf(tmp_path / "t.vcf.gz", reference, _sorted(records, reference))
    sim = SimulationConfig(profile="p.yaml", flank_bp=kw.pop("flank_bp", 400))  # type: ignore[arg-type]
    return prepare_sample(
        sample_id="s1",
        sex=Sex.MALE,
        truth_vcf=truth,
        planted_ids=planted_ids(records),
        reference=reference,
        pars=PARS,
        sim=sim,
        extra_bed=None,
        min_segment_bp=300,
        sim_dir=tmp_path / "sim",
        **kw,  # type: ignore[arg-type]
    )


def _mason_records(path: Path) -> list[tuple[str, int, str, str, str]]:
    out = []
    with pysam.VariantFile(str(path)) as vcf:
        for r in vcf:
            gt = r.samples["s1"]["GT"]
            sep = "|" if r.samples["s1"].phased else "/"
            out.append((r.chrom, r.pos, r.ref, r.alts[0], sep.join(str(a) for a in gt)))
    return out


def test_prepare_sample_writes_valid_mason_inputs(
    tmp_path: Path, reference: Reference, records: list[Rec]
) -> None:
    """GT only, single ploidy per file, contigs = FASTA, decomposed, shifted, overlaps out."""
    prep = _prepare(tmp_path, reference, records)
    assert [g.group for g in prep.groups] == ["diploid", "haploid"]
    assert prep.decomposed == [planted_ids(records)[0]]
    reasons = sorted((r.pos, reason) for r, reason in prep.dropped)
    assert reasons == [(1002, "overlap_on_haplotype"), (3100, "het_on_haploid_segment")]
    diploid, haploid = prep.groups
    assert [s.name for s in diploid.segments] == ["chr1_600_1402", "chr2_3600_4400"]
    assert [s.name for s in haploid.segments] == ["chrX_2600_3400"]
    assert diploid.vcf is not None and haploid.vcf is not None
    dip = _mason_records(diploid.vcf)
    # chr1:1000 delins (3 bp -> 2 bp, different first base) -> SNV at 1000 + SNV 1001 + deletion
    # anchored at 1001; background SNV at 1001 on haplotype 1 is kept (other haplotype),
    # background SNV at 1002 on haplotype 2 overlaps the planted deletion and is dropped.
    assert [(c, p, g) for c, p, _, _, g in dip] == [
        ("chr1_600_1402", 401, "0|1"),
        ("chr1_600_1402", 402, "1|0"),
        ("chr1_600_1402", 402, "0|1"),
        ("chr1_600_1402", 402, "0|1"),
        ("chr1_600_1402", 701, "1|1"),
        ("chr2_3600_4400", 401, "1|0"),
    ]
    hap = _mason_records(haploid.vcf)
    assert [(c, p, g) for c, p, _, _, g in hap] == [
        ("chrX_2600_3400", 401, "1"),
        ("chrX_2600_3400", 601, "1"),  # 1|1 background on a haploid segment -> 1
    ]
    assert validate_mason_vcf(diploid.vcf, diploid.fasta) == 6
    with pysam.FastaFile(str(diploid.fasta)) as fa:
        assert list(fa.references) == ["chr1_600_1402", "chr2_3600_4400"]
        assert fa.fetch("chr1_600_1402") == reference.fetch("chr1", 599, 1402)
    with pysam.VariantFile(str(diploid.vcf)) as vcf:
        assert list(vcf.header.formats) == ["GT"]
    assert prep.segments_tsv.read_text().splitlines()[1].startswith("chr1_600_1402\tchr1\t600")
    assert [r.variant_id for r in prep.planted] == planted_ids(records)


def test_background_delins_with_different_first_base(tmp_path: Path, reference: Reference) -> None:
    """A delins whose first base differs decomposes into an SNV and an indel at the same POS;
    the SNV must be written first or the Mason VCF fails validation (HG002, e.g. chr17:29305771
    C>TGGAGATGGGGTTT)."""
    r = reference
    records: list[Rec] = [
        ("chr1", 1000, _ref(r, "chr1", 1000), _other(_ref(r, "chr1", 1000)) + "GG", (1, 1), None),
        ("chr1", 1100, _ref(r, "chr1", 1100), _other(_ref(r, "chr1", 1100)), (1, 0), "p:snv"),
        ("chr1", 1200, _ref(r, "chr1", 1200, 3), _other(_ref(r, "chr1", 1200)), (0, 1), None),
    ]
    prep = _prepare(tmp_path, reference, records)
    assert prep.dropped == []
    (diploid,) = prep.groups
    assert diploid.vcf is not None
    dip = _mason_records(diploid.vcf)
    kinds = [(p, "snv" if len(ref) == len(alt) else "indel") for _, p, ref, alt, _ in dip]
    assert kinds == [(301, "snv"), (301, "indel"), (401, "snv"), (501, "snv"), (501, "indel")]
    assert validate_mason_vcf(diploid.vcf, diploid.fasta) == 5


def test_planted_overlap_is_an_error(tmp_path: Path, reference: Reference) -> None:
    r = reference
    records: list[Rec] = [
        ("chr1", 1000, _ref(r, "chr1", 1000, 4), _ref(r, "chr1", 1000), (1, 0), "p:a"),
        ("chr1", 1002, _ref(r, "chr1", 1002), _other(_ref(r, "chr1", 1002)), (1, 0), "p:b"),
    ]
    with pytest.raises(SimPrepError, match="overlaps another planted variant"):
        _prepare(tmp_path, reference, records)


def test_short_segments_are_an_error(tmp_path: Path, reference: Reference) -> None:
    r = reference
    records: list[Rec] = [
        (
            "chr1",
            1000,
            _ref(r, "chr1", 1000),
            "T" if _ref(r, "chr1", 1000) != "T" else "G",
            (1, 0),
            "p",
        )
    ]
    with pytest.raises(SimPrepError, match="shorter than the longest fragment"):
        _prepare(tmp_path, reference, records, flank_bp=50)


def test_genome_mode_includes_all_background(
    tmp_path: Path, reference: Reference, records: list[Rec]
) -> None:
    truth = _truth_vcf(tmp_path / "t.vcf.gz", reference, _sorted(records, reference))
    prep = prepare_sample(
        sample_id="s1",
        sex=Sex.MALE,
        truth_vcf=truth,
        planted_ids=planted_ids(records),
        reference=reference,
        pars=PARS,
        sim=SimulationConfig(profile="p.yaml", mode=SimulationMode.GENOME),
        extra_bed=None,
        min_segment_bp=300,
        sim_dir=tmp_path / "sim",
    )
    names = [s.name for g in prep.groups for s in g.segments]
    assert names == [
        "chr1_1_10000",
        "chr2_1_10000",
        "chrX_1_1000",
        "chrX_1001_8000",
        "chrY_501_2000",
    ]
    assert len(read_truth(truth)) == len(records)


def _seg(name: str, length: int) -> SegmentRecord:
    return SegmentRecord(name=name, chrom="chr1", start=1, end=length, ploidy=2)


def test_batch_segments_packs_consecutive_segments() -> None:
    segs = [_seg("a", 60), _seg("b", 30), _seg("c", 20), _seg("d", 150), _seg("e", 10)]
    assert batch_segments(segs, None) == [segs]
    names = [[s.name for s in b] for b in batch_segments(segs, 100)]
    assert names == [["a", "b"], ["c"], ["d"], ["e"]]  # d is longer than the limit: alone
    assert [[s.name for s in b] for b in batch_segments(segs, 1000)] == [["a", "b", "c", "d", "e"]]


def test_genome_mode_batches(tmp_path: Path, reference: Reference, records: list[Rec]) -> None:
    truth = _truth_vcf(tmp_path / "t.vcf.gz", reference, _sorted(records, reference))

    def prep(sim_dir: str, batch_bp: int | None):  # type: ignore[no-untyped-def]
        return prepare_sample(
            sample_id="s1",
            sex=Sex.MALE,
            truth_vcf=truth,
            planted_ids=planted_ids(records),
            reference=reference,
            pars=PARS,
            sim=SimulationConfig(
                profile="p.yaml", mode=SimulationMode.GENOME, mason_batch_bp=batch_bp
            ),
            extra_bed=None,
            min_segment_bp=300,
            sim_dir=tmp_path / sim_dir,
        )

    whole = prep("whole", None)
    batched = prep("batched", 15_000)
    assert [(g.name, g.batched) for g in whole.groups] == [("diploid", False), ("haploid", False)]
    assert [(g.name, [s.name for s in g.segments]) for g in batched.groups] == [
        ("diploid-1", ["chr1_1_10000"]),
        ("diploid-2", ["chr2_1_10000", "chrX_1_1000"]),
        ("haploid-1", ["chrX_1001_8000", "chrY_501_2000"]),
    ]
    for g in batched.groups:
        assert g.batched and g.fasta.name == f"{g.name}.fa"
        assert g.vcf is not None and g.vcf.name == f"{g.name}.mason.vcf"
    # The batches carry exactly the records of the unbatched groups.
    by_group: dict[str, list[tuple[str, int, str, str, str]]] = {}
    for g in batched.groups:
        assert g.vcf is not None
        by_group.setdefault(g.group, []).extend(_mason_records(g.vcf))
    for g in whole.groups:
        assert g.vcf is not None
        assert by_group[g.group] == _mason_records(g.vcf)
        assert sum(b.records for b in batched.groups if b.group == g.group) == g.records


def _write_mason(
    path: Path,
    segments: list[SegmentRecord],
    rows: list[tuple[str, int, str, str, tuple[int, ...]]],
) -> Path:
    writer = MasonVcfWriter(path, "s1", segments)
    for seg, pos, ref, alt, gt in rows:
        writer.write(seg, pos, ref, alt, gt, ".")
    writer.close()
    return path


def test_validate_mason_vcf_rejections(tmp_path: Path) -> None:
    fasta = tmp_path / "g.fa"
    fasta.write_text(">s1\n" + "ACGT" * 50 + "\n>s2\n" + "ACGT" * 50 + "\n")
    pysam.faidx(str(fasta))
    segs = [
        SegmentRecord(name="s1", chrom="chr1", start=1, end=200, ploidy=2),
        SegmentRecord(name="s2", chrom="chr1", start=301, end=500, ploidy=2),
    ]
    good = _write_mason(tmp_path / "good.vcf", segs, [("s1", 10, "G", "T", (0, 1))])
    assert validate_mason_vcf(good, fasta) == 1
    cases = {
        "MNV/delins": [("s1", 10, "GT", "AA", (0, 1))],
        "ploidy": [("s1", 10, "G", "T", (0, 1)), ("s1", 20, "A", "T", (1,))],
        "conflict": [("s1", 10, "GTA", "G", (1, 0)), ("s1", 11, "T", "C", (1, 0))],
        "not sorted": [("s2", 10, "G", "T", (0, 1)), ("s1", 10, "G", "T", (0, 1))],
    }
    for message, rows in cases.items():
        path = _write_mason(tmp_path / f"{message[:3]}.vcf", segs, rows)
        with pytest.raises(SimPrepError, match=message):
            validate_mason_vcf(path, fasta)
    reversed_segs = list(reversed(segs))
    path = _write_mason(tmp_path / "order.vcf", reversed_segs, [])
    with pytest.raises(SimPrepError, match="header contigs differ"):
        validate_mason_vcf(path, fasta)
    ps = tmp_path / "ps.vcf"
    text = good.read_text().replace("\tGT\t0|1", "\tGT:PS\t0|1:10")
    text = text.replace(
        "##FORMAT=<ID=GT", '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="p">\n##FORMAT=<ID=GT'
    )
    ps.write_text(text)
    with pytest.raises(SimPrepError, match="only GT"):
        validate_mason_vcf(ps, fasta)


def test_boundary_crossing_background_logged_once(tmp_path: Path, reference: Reference) -> None:
    r = reference
    snv = ("chr1", 1000, _ref(r, "chr1", 1000), _other(_ref(r, "chr1", 1000)), (1, 0), "p")
    crossing = ("chr1", 595, _ref(r, "chr1", 595, 12), _ref(r, "chr1", 595), (0, 1), None)
    prep = _prepare(tmp_path, reference, [crossing, snv])
    assert [(rec.pos, reason) for rec, reason in prep.dropped] == [
        (595, "crosses_segment_boundary")
    ]
    assert prep.dropped_log is not None
    assert prep.dropped_log.read_text().splitlines()[1].endswith("crosses_segment_boundary")


def test_indel_at_segment_end_is_dropped(tmp_path: Path, reference: Reference) -> None:
    """mason_simulator crashes on a deletion that reaches the last base of a segment."""
    r = reference
    snv = ("chr1", 1000, _ref(r, "chr1", 1000), _other(_ref(r, "chr1", 1000)), (1, 0), "p")
    deletion = ("chr1", 1398, _ref(r, "chr1", 1398, 3), _ref(r, "chr1", 1398), (1, 1), None)
    last_snv = ("chr1", 1400, _ref(r, "chr1", 1400), _other(_ref(r, "chr1", 1400)), (0, 1), None)
    prep = _prepare(tmp_path, reference, [snv, deletion, last_snv])
    assert [s.end for s in prep.segments] == [1400]
    assert [(rec.pos, reason) for rec, reason in prep.dropped] == [(1398, "indel_at_segment_end")]


def test_missing_planted_id_is_an_error(tmp_path: Path, reference: Reference) -> None:
    r = reference
    records: list[Rec] = [("chr1", 1000, _ref(r, "chr1", 1000), "T", (1, 0), "p")]
    truth = _truth_vcf(tmp_path / "t.vcf.gz", reference, records)
    with pytest.raises(SimPrepError, match="planted variants not found: curated:chr1:999"):
        prepare_sample(
            sample_id="s1",
            sex=Sex.MALE,
            truth_vcf=truth,
            planted_ids=["curated:chr1:999:A:T"],
            reference=reference,
            pars=PARS,
            sim=SimulationConfig(profile="p.yaml"),
            extra_bed=None,
            min_segment_bp=300,
            sim_dir=tmp_path / "sim",
        )
