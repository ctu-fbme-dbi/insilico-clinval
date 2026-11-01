"""Variant-type classification and allele sanity checks."""

from __future__ import annotations

import re
from typing import Protocol

from insilico_clinval.schemas.enums import VariantType

_BASES = re.compile(r"^[ACGT]+$")

INDEL_TYPES = (VariantType.INSERTION, VariantType.DELETION, VariantType.DELINS)


def is_simple_allele(allele: str | None) -> bool:
    """True for a plain nucleotide allele (no symbolic, missing, spanning or N alleles)."""
    return allele is not None and bool(_BASES.match(allele.upper()))


def classify(ref: str, alt: str) -> VariantType:
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        raise ValueError(f"REF equals ALT: {ref}")
    if len(ref) == len(alt):
        return VariantType.SNV if len(ref) == 1 else VariantType.MNV
    if len(ref) == 1 and alt[0] == ref:
        return VariantType.INSERTION
    if len(alt) == 1 and ref[0] == alt:
        return VariantType.DELETION
    return VariantType.DELINS


def indel_length(ref: str, alt: str) -> int:
    """``len(alt) - len(ref)`` (docs/schemas.md); 0 for SNV/MNV."""
    return len(alt) - len(ref)


def event_length(ref: str, alt: str) -> int:
    """Number of inserted/deleted/replaced bases, used by the indel ladder and max length.

    Insertion: ``len(alt) - 1``; deletion: ``len(ref) - 1``; delins: longer allele minus the
    shared anchor base (if any); SNV 1, MNV ``len(ref)``.
    """
    vtype = classify(ref, alt)
    if vtype in (VariantType.SNV, VariantType.MNV):
        return len(ref)
    anchor = 1 if ref[0].upper() == alt[0].upper() else 0
    return max(len(ref), len(alt)) - anchor


def ref_span0(pos: int, ref: str) -> tuple[int, int]:
    """Half-open 0-based reference span of REF."""
    return pos - 1, pos - 1 + len(ref)


class SequenceSource(Protocol):
    def fetch(self, chrom: str, start0: int, end0: int) -> str: ...


def left_normalize(
    chrom: str, pos: int, ref: str, alt: str, reference: SequenceSource
) -> tuple[int, str, str]:
    """Left-align and trim a biallelic variant (Tan et al. 2015, Bioinformatics 31:2202).

    Repeatedly truncate a shared last base and, when an allele becomes empty, extend both alleles
    by the preceding reference base; then trim shared leading bases while both alleles are longer
    than one base. The result is what ``bcftools norm -f`` produces for a biallelic record.
    """
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        raise ValueError(f"{chrom}:{pos}: REF equals ALT ({ref})")
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif not ref or not alt:
            if pos == 1:
                raise ValueError(f"{chrom}:{pos}: cannot anchor an indel at the contig start")
            pos -= 1
            base = reference.fetch(chrom, pos - 1, pos)
            ref, alt = base + ref, base + alt
        else:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt
