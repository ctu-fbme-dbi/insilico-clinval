"""Reference-build comparison: primary contigs and GIAB-style masked regions."""

from __future__ import annotations

from pathlib import Path

import pysam
from typer.testing import CliRunner

from insilico_clinval.cli import app
from insilico_clinval.genome import refdiff
from insilico_clinval.genome.reference import Reference


def _fasta(path: Path, contigs: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in contigs))
    pysam.faidx(str(path))
    return path


UNMASKED = [("chr1", "ACGT" * 30), ("chr2", "TTGCA" * 19 + "TTRCA")]  # IUPAC R at chr2:98
# chr1 bases 10-19 and 115-119 masked; chr2: only the IUPAC code turned into N (not counted)
MASKED = [
    ("chr1", "ACGTACGTAC" + "N" * 10 + ("ACGT" * 30)[20:115] + "N" * 5),
    ("chr2", "TTGCA" * 19 + "TTNCA"),
]


def test_primary_problems(tmp_path: Path) -> None:
    a = Reference(_fasta(tmp_path / "a.fa", UNMASKED))
    same = Reference(_fasta(tmp_path / "b.fa", [*MASKED, ("chrUn_decoy", "ACGT")]))
    assert refdiff.primary_problems(a, same) == []  # extra decoy contigs are ignored
    short = Reference(_fasta(tmp_path / "c.fa", [("chr1", "ACGT" * 29), UNMASKED[1]]))
    assert refdiff.primary_problems(a, short) == ["chr1: length 120 vs 116"]
    swapped = Reference(_fasta(tmp_path / "d.fa", [UNMASKED[1], UNMASKED[0]]))
    assert refdiff.primary_problems(a, swapped) == ["primary contigs are in a different order"]
    missing = Reference(_fasta(tmp_path / "e.fa", [UNMASKED[0]]))
    assert refdiff.primary_problems(a, missing)[0].startswith("chr2: missing from")


def test_masked_regions(tmp_path: Path) -> None:
    masked = Reference(_fasta(tmp_path / "m.fa", MASKED))
    plain = Reference(_fasta(tmp_path / "u.fa", UNMASKED))
    assert list(refdiff.masked_only_in(masked, plain)) == [("chr1", 10, 20), ("chr1", 115, 120)]
    # runs spanning a chunk border are merged
    refdiff_chunk = refdiff.CHUNK_BP
    try:
        refdiff.CHUNK_BP = 12
        assert list(refdiff.masked_only_in(masked, plain)) == [
            ("chr1", 10, 20), ("chr1", 115, 120)
        ]  # fmt: skip
    finally:
        refdiff.CHUNK_BP = refdiff_chunk
    assert list(refdiff.masked_only_in(plain, masked)) == []


def test_check_references_cli(tmp_path: Path) -> None:
    m, u = _fasta(tmp_path / "m.fa", MASKED), _fasta(tmp_path / "u.fa", UNMASKED)
    bed = tmp_path / "masked.bed"
    result = CliRunner().invoke(
        app,
        ["check-references", "--reference", str(m), "--against", str(u), "--masked-bed", str(bed)],
    )
    assert result.exit_code == 0, result.output
    assert bed.read_text() == "chr1\t10\t20\nchr1\t115\t120\n"
    bad = _fasta(tmp_path / "s.fa", [("chr1", "ACGT"), UNMASKED[1]])
    result = CliRunner().invoke(
        app, ["check-references", "--reference", str(m), "--against", str(bad)]
    )
    assert result.exit_code == 1
