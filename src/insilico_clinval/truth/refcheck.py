"""REF alleles against the reference FASTA."""

from __future__ import annotations

from collections.abc import Iterable

from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.variant import TestVariant
from insilico_clinval.truth import TruthError


def check_ref_alleles(variants: Iterable[TestVariant], reference: Reference) -> None:
    """Raise TruthError listing every variant whose REF does not match the FASTA."""
    errors: list[str] = []
    for v in variants:
        length = reference.lengths.get(v.chrom)
        if length is None:
            errors.append(f"{v.variant_id}: contig {v.chrom} not in reference")
            continue
        end = v.pos - 1 + len(v.ref)
        if end > length:
            errors.append(f"{v.variant_id}: REF extends beyond {v.chrom} (length {length})")
            continue
        actual = reference.fetch(v.chrom, v.pos - 1, end)
        if actual != v.ref.upper():
            errors.append(
                f"{v.variant_id}: REF {v.ref} does not match reference {actual} "
                f"at {v.chrom}:{v.pos} (source {v.source.database} {v.source.record_id})"
            )
    if errors:
        raise TruthError(
            "REF allele mismatch against the reference FASTA:\n  " + "\n  ".join(errors)
        )
