"""Minimal genomic HGVS (``NC_…:g.``) parsing and conversion to VCF alleles.

Only the forms needed for expert-curated SNVs and small indels are supported: substitution
(``g.100A>G``), deletion (``g.100del``, ``g.100_102delACT``), duplication (``g.100_102dup``),
insertion (``g.100_101insAC``) and deletion-insertion (``g.100_102delinsTT``). Uncertain
positions (``(…)``, ``?``), intronic offsets, repeats (``[n]``), inversions, conversions and
identity (``=``) raise :class:`HgvsUnsupported`; callers skip such records and count them.

Reference bases stated in the HGVS (substitution REF, optional deleted / duplicated sequence) are
checked against the FASTA; a difference raises :class:`HgvsMismatch`. The returned VCF alleles are
anchored as in VCF but *not* normalised (HGVS is 3'-shifted); use
:func:`insilico_clinval.genome.alleles.left_normalize`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal

from insilico_clinval.genome.alleles import SequenceSource

HgvsKind = Literal["sub", "del", "dup", "ins", "delins"]

_GENOMIC = re.compile(r"^(?P<acc>N[CTW]_\d+\.\d+):g\.(?P<start>\d+)(?:_(?P<end>\d+))?(?P<edit>.+)$")
_EDITS: tuple[tuple[HgvsKind, re.Pattern[str]], ...] = (
    ("sub", re.compile(r"^(?P<ref>[ACGT])>(?P<alt>[ACGT])$")),
    ("delins", re.compile(r"^delins(?P<alt>[ACGT]+)$")),
    ("del", re.compile(r"^del(?P<ref>[ACGT]*)$")),
    ("dup", re.compile(r"^dup(?P<ref>[ACGT]*)$")),
    ("ins", re.compile(r"^ins(?P<alt>[ACGT]+)$")),
)


class HgvsUnsupported(ValueError):
    """The expression is valid HGVS (or not HGVS at all) but outside the supported subset."""


class HgvsMismatch(ValueError):
    """Reference bases stated in the expression differ from the FASTA."""


@dataclass(frozen=True)
class GenomicHgvs:
    accession: str
    kind: HgvsKind
    start: int
    end: int
    ref: str | None  # stated reference bases (substitution REF, del/dup sequence), if any
    alt: str | None  # substitution ALT, inserted sequence

    @property
    def text(self) -> str:
        span = str(self.start) if self.start == self.end else f"{self.start}_{self.end}"
        if self.kind == "sub":
            return f"{self.accession}:g.{span}{self.ref}>{self.alt}"
        stated = (self.ref or "") if self.kind in ("del", "dup") else (self.alt or "")
        return f"{self.accession}:g.{span}{self.kind}{stated}"


def parse_genomic(expression: str) -> GenomicHgvs:
    text = expression.strip()
    m = _GENOMIC.match(text)
    if m is None:
        raise HgvsUnsupported(f"not a supported genomic HGVS expression: {text!r}")
    start = int(m["start"])
    end = int(m["end"]) if m["end"] else start
    if end < start:
        raise HgvsUnsupported(f"{text}: end before start")
    edit = m["edit"]
    for kind, pattern in _EDITS:
        e = pattern.match(edit)
        if e is None:
            continue
        groups = e.groupdict()
        ref = groups.get("ref") or None
        alt = groups.get("alt") or None
        if kind == "sub" and end != start:
            raise HgvsUnsupported(f"{text}: substitution over a range")
        if kind == "ins" and end != start + 1:
            raise HgvsUnsupported(f"{text}: insertion must be between adjacent positions")
        if kind in ("del", "dup") and ref is not None and len(ref) != end - start + 1:
            raise HgvsUnsupported(f"{text}: stated sequence length differs from the range")
        return GenomicHgvs(m["acc"], kind, start, end, ref, alt)
    raise HgvsUnsupported(f"{text}: unsupported edit {edit!r}")


def to_vcf(
    h: GenomicHgvs, chrom: str, reference: SequenceSource, length: int
) -> tuple[int, str, str]:
    """VCF ``(pos, ref, alt)`` of the expression on ``chrom`` (anchored, not normalised)."""
    if h.end > length:
        raise HgvsUnsupported(f"{h.text}: beyond the end of {chrom} (length {length})")
    span = reference.fetch(chrom, h.start - 1, h.end)
    if h.ref is not None and span != h.ref:
        raise HgvsMismatch(
            f"{h.text}: stated reference {h.ref} differs from the FASTA ({span}) at "
            f"{chrom}:{h.start}"
        )
    if h.kind == "sub":
        assert h.alt is not None
        return h.start, span, h.alt
    if h.kind == "delins":
        assert h.alt is not None
        if span == h.alt:
            raise HgvsUnsupported(f"{h.text}: inserted sequence equals the deleted one")
        return h.start, span, h.alt
    if h.kind == "del":
        if h.start == 1:
            raise HgvsUnsupported(f"{h.text}: deletion at the contig start")
        anchor = reference.fetch(chrom, h.start - 2, h.start - 1)
        return h.start - 1, anchor + span, anchor
    if h.kind == "dup":
        anchor = span[-1]
        return h.end, anchor, anchor + span
    assert h.alt is not None  # ins
    anchor = reference.fetch(chrom, h.start - 1, h.start)
    return h.start, anchor, anchor + h.alt
