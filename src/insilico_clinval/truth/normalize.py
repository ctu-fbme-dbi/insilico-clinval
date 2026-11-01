"""Normalisation with ``bcftools norm`` and the idempotence assertion.

Selection works on source coordinates. The truth VCF is passed through
``bcftools norm -f REF -c e -m -any``; if any planted record changes (i.e. the selected
representation was not left-aligned / biallelic), deduplication and variant ids would be stale,
so this is an error. Background records (INFO ``TV_BG``) may be split and realigned; the output
must remain coordinate-sorted in FASTA-index order.
"""

from __future__ import annotations

from pathlib import Path

import pysam

from insilico_clinval._runner import Toolbox
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.truth import TruthError


def _planted(path: Path) -> list[tuple[str, int, str | None, tuple[str, ...]]]:
    """Planted records (not background) of a VCF, streamed; also asserts coordinate order."""
    out: list[tuple[str, int, str | None, tuple[str, ...]]] = []
    with pysam.VariantFile(str(path)) as vcf:
        order = {str(c): i for i, c in enumerate(vcf.header.contigs)}
        previous = (-1, 0)
        for r in vcf:
            key = (order[r.chrom], r.pos)
            if key < previous:
                raise TruthError(f"{path.name}: records out of order at {r.chrom}:{r.pos}")
            previous = key
            if "TV_BG" not in r.info:
                out.append((r.chrom, r.pos, r.ref, tuple(r.alts or ())))
    return out


def normalize_vcf(
    src: Path,
    dst: Path,
    fasta: Path,
    *,
    step: str,
    manifest: ManifestBuilder | None,
    toolbox: Toolbox,
) -> Path:
    """Run bcftools norm on ``src`` into ``dst`` (plain VCF) and assert it was a no-op."""
    args = ["norm", "--no-version", "-f", str(fasta.resolve()), "-c", "e", "-m", "-any"]
    toolbox.run(
        "bcftools",
        [*args, "-O", "v", "-o", str(dst.resolve()), str(src.resolve())],
        step=step,
        manifest=manifest,
        mounts=[src.parent, dst.parent, fasta.parent],
    )
    before, after = _planted(src), _planted(dst)
    if before != after:
        changed = [f"{b} -> {a}" for b, a in zip(before, after, strict=False) if b != a]
        if len(before) != len(after):
            changed.append(f"record count {len(before)} -> {len(after)}")
        raise TruthError(
            "bcftools norm changed the truth set (selected variants are not normalised; "
            "fix the source or the explicit variant):\n  " + "\n  ".join(changed)
        )
    return dst


def compress_and_index(src: Path, dst_gz: Path) -> Path:
    """bgzip + tabix via pysam (deterministic)."""
    pysam.tabix_compress(str(src), str(dst_gz), force=True)
    pysam.tabix_index(str(dst_gz), preset="vcf", force=True)
    return dst_gz
