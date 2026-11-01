"""Comparison of two builds of the same reference (e.g. GRCh38 no-alt vs GIAB-masked GRCh38).

A suite may simulate from one FASTA (the suite reference) while the pipeline aligns to another
(the adapter's reference). Coordinates are shared only if the primary contigs agree in name,
order and length; extra contigs (decoys, ALT, HLA) are ignored. Positions that are ``N`` in one
FASTA and a plain base (ACGT) in the other (e.g. the false duplications of GRCh38 masked by GIAB)
form a BED; IUPAC ambiguity codes turned into ``N`` are not counted.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from insilico_clinval.genome.reference import Reference

PRIMARY_CONTIGS: tuple[str, ...] = (*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY", "chrM")

CHUNK_BP = 1_000_000


def primary_problems(a: Reference, b: Reference) -> list[str]:
    """Differences of the primary contigs (missing, length, relative order); empty if equal."""
    problems: list[str] = []
    for name in PRIMARY_CONTIGS:
        la, lb = a.lengths.get(name), b.lengths.get(name)
        if la is None and lb is None:
            continue
        if la is None or lb is None:
            missing = a.fasta if la is None else b.fasta
            problems.append(f"{name}: missing from {missing}")
        elif la != lb:
            problems.append(f"{name}: length {la} vs {lb}")
    shared = [n for n in PRIMARY_CONTIGS if n in a.lengths and n in b.lengths]
    if sorted(shared, key=a.order.__getitem__) != sorted(shared, key=b.order.__getitem__):
        problems.append("primary contigs are in a different order")
    return problems


def masked_only_in(
    masked: Reference, other: Reference, contigs: tuple[str, ...] = PRIMARY_CONTIGS
) -> Iterator[tuple[str, int, int]]:
    """Half-open intervals that are ``N`` in ``masked`` and ACGT in ``other`` (same lengths)."""
    for name in contigs:
        if masked.lengths.get(name) != other.lengths.get(name) or name not in masked.lengths:
            continue
        length = masked.lengths[name]
        run_start: int | None = None
        for start in range(0, length, CHUNK_BP):
            end = min(start + CHUNK_BP, length)
            a, b = masked.fetch(name, start, end), other.fetch(name, start, end)
            if "N" not in a or a == b:
                if run_start is not None:
                    yield name, run_start, start
                    run_start = None
                continue
            for i, (x, y) in enumerate(zip(a, b, strict=True)):
                hit = x == "N" and y in "ACGT"
                if hit and run_start is None:
                    run_start = start + i
                elif not hit and run_start is not None:
                    yield name, run_start, start + i
                    run_start = None
        if run_start is not None:
            yield name, run_start, length


def write_masked_bed(masked: Reference, other: Reference, path: Path) -> int:
    """Write the BED of :func:`masked_only_in`; returns the number of bases."""
    total = 0
    with path.open("w", encoding="utf-8") as fh:
        for chrom, start, end in masked_only_in(masked, other):
            fh.write(f"{chrom}\t{start}\t{end}\n")
            total += end - start
    return total
