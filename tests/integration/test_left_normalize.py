"""left_normalize agrees with bcftools norm on random indels of the fixture reference."""

from __future__ import annotations

import random
from pathlib import Path

import pysam
import pytest

from insilico_clinval._runner import Toolbox
from insilico_clinval.genome.alleles import left_normalize
from insilico_clinval.genome.reference import Reference
from tests.conftest import FIXTURES, requires_bcftools

pytestmark = [pytest.mark.integration, requires_bcftools]


def _random_events(reference: Reference, n: int) -> list[tuple[str, int, str, str]]:
    rng = random.Random(7)
    events: list[tuple[str, int, str, str]] = []
    for i in range(n):
        chrom = "chr1" if i % 2 else "chr2"
        pos = rng.randrange(50, reference.lengths[chrom] - 50)
        length = rng.randint(1, 12)
        span = reference.fetch(chrom, pos - 1, pos + length)
        kind = i % 3
        if kind == 0:  # deletion, possibly 3'-shifted inside a repeat
            ref, alt = span, span[0]
        elif kind == 1:  # duplication of the following bases (tandem repeat)
            ref, alt = span[0], span[0] + span[1:]
        else:  # delins with a shared first base
            ref, alt = (
                span[: length // 2 + 2],
                span[0] + "".join(rng.choice("ACGT") for _ in range(3)),
            )
            if ref == alt:
                continue
        npos, nref, nalt = left_normalize(chrom, pos, ref, alt, reference)
        events.append((chrom, npos, nref, nalt))
    return events


def test_matches_bcftools(tmp_path: Path, reference: Reference) -> None:
    events = sorted(
        set(_random_events(reference, 300)), key=lambda e: (reference.order[e[0]], e[1], e[2], e[3])
    )
    header = pysam.VariantHeader()
    for name, length in reference.contigs:
        header.contigs.add(name, length=length)
    src = tmp_path / "in.vcf"
    with pysam.VariantFile(str(src), "w", header=header) as out:
        for chrom, pos, ref, alt in events:
            out.write(out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt)))
    dst = tmp_path / "out.vcf"
    fasta = FIXTURES / "ref.fa"
    Toolbox().run(
        "bcftools",
        [
            "norm",
            "--no-version",
            "-f",
            str(fasta),
            "-c",
            "e",
            "-O",
            "v",
            "-o",
            str(dst),
            str(src),
        ],
        step="test",
        manifest=None,
    )
    with pysam.VariantFile(str(dst)) as vcf:
        normalised = sorted(
            ((r.chrom, r.pos, str(r.ref), str((r.alts or ("",))[0])) for r in vcf),
            key=lambda e: (reference.order[e[0]], e[1], e[2], e[3]),
        )
    assert normalised == events
