"""Genomic HGVS -> VCF conversion and left normalisation (eRepo coordinates)."""

from __future__ import annotations

import pytest

from insilico_clinval.genome.alleles import left_normalize
from insilico_clinval.genome.hgvs import HgvsMismatch, HgvsUnsupported, parse_genomic, to_vcf
from insilico_clinval.genome.reference import Reference

# Fixture bases (tests/fixtures/README.md): chr2:9499 C, chr2:9500-9505 AAAAAA, chr2:9506 G,
# chr2:5100 G, chr2:9800 C.
CHR2 = "NC_000002.12"


def _vcf(expr: str, reference: Reference) -> tuple[int, str, str]:
    h = parse_genomic(expr)
    pos, ref, alt = to_vcf(h, "chr2", reference, reference.lengths["chr2"])
    return left_normalize("chr2", pos, ref, alt, reference)


@pytest.mark.parametrize(
    ("expr", "expected"),
    [
        (f"{CHR2}:g.5100G>T", (5100, "G", "T")),
        # HGVS is 3'-shifted; VCF is left-aligned before the homopolymer chr2:9500-9505
        (f"{CHR2}:g.9505del", (9499, "CA", "C")),
        (f"{CHR2}:g.9505delA", (9499, "CA", "C")),
        (f"{CHR2}:g.9504_9505del", (9499, "CAA", "C")),
        (f"{CHR2}:g.9505dup", (9499, "C", "CA")),
        (f"{CHR2}:g.9502_9503insA", (9499, "C", "CA")),
        (f"{CHR2}:g.9500_9505delAAAAAA", (9499, "CAAAAAA", "C")),
        (f"{CHR2}:g.9505_9506delinsTT", (9505, "AG", "TT")),
        (f"{CHR2}:g.9506delinsAT", (9506, "G", "AT")),
    ],
)
def test_hgvs_to_normalised_vcf(
    expr: str, expected: tuple[int, str, str], reference: Reference
) -> None:
    assert _vcf(expr, reference) == expected


def test_stated_reference_is_checked(reference: Reference) -> None:
    with pytest.raises(HgvsMismatch, match=r"stated reference A differs from the FASTA \(C\)"):
        _vcf(f"{CHR2}:g.9800A>G", reference)
    with pytest.raises(HgvsMismatch):
        _vcf(f"{CHR2}:g.9505delC", reference)


@pytest.mark.parametrize(
    "expr",
    [
        f"{CHR2}:g.(9400_9450)_(9500_9550)del",
        f"{CHR2}:g.9500_9505inv",
        f"{CHR2}:g.9500A[8]",
        f"{CHR2}:g.9500=",
        f"{CHR2}:g.9500_9505ins",
        f"{CHR2}:g.9500_9502insA",
        f"{CHR2}:g.9500_9501delA",
        "NM_000001.1:c.100A>G",
        "not hgvs",
    ],
)
def test_unsupported_expressions(expr: str, reference: Reference) -> None:
    with pytest.raises(HgvsUnsupported):
        _vcf(expr, reference)


def test_left_normalize_trims_and_shifts(reference: Reference) -> None:
    # MNV-like representation with shared flanks collapses to the minimal SNV
    seq = reference.fetch("chr2", 5098, 5101)  # chr2:5099-5101
    pos, ref, alt = left_normalize("chr2", 5099, seq, seq[0] + "T" + seq[2], reference)
    assert (pos, ref, alt) == (5100, "G", "T")
    with pytest.raises(ValueError, match="REF equals ALT"):
        left_normalize("chr2", 5100, "G", "G", reference)
