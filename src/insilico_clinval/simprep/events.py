"""Decomposition of alleles into Mason-representable parts and per-haplotype overlap checks.

``mason_simulator`` (seqan-v2.4.0, ``vcf_materialization.cpp``) materialises three kinds of
small records: SNVs (``len(REF) == len(ALT) == 1``), insertions (``len(REF) == 1``; the bases after
the first ALT base are inserted after POS) and deletions (``len(ALT) == 1``; ``len(REF) - 1``
bases after POS are deleted). Records with both alleles longer than one base (MNV, delins) are
skipped silently, so they are decomposed here: SNVs for the differing bases of the common
prefix, plus one insertion or deletion anchored at the last common base. Applying the parts
reproduces the ALT haplotype exactly (the anchor base of an indel part is never read by Mason).

Mason applies the parts of one haplotype in position order without any overlap check, so every
haplotype must be free of conflicts (:func:`conflicts`).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class Part:
    """A Mason-representable record: 1-based POS, REF, ALT (SNV, insertion or deletion)."""

    pos: int
    ref: str
    alt: str

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)

    def span0(self) -> tuple[int, int]:
        """Half-open 0-based reference span changed by the part (empty for an insertion)."""
        if not self.is_indel:
            return self.pos - 1, self.pos
        # Indels are anchored at POS; the change starts after the anchor base.
        return self.pos, self.pos - 1 + len(self.ref)


def decompose(pos: int, ref: str, alt: str) -> list[Part]:
    """Split one REF/ALT allele into Mason-representable parts (sorted)."""
    ref, alt = ref.upper(), alt.upper()
    if ref == alt:
        raise ValueError(f"REF equals ALT at {pos}: {ref}")
    m, n = len(ref), len(alt)
    if (m == 1 or n == 1) and ref[0] == alt[0] and m != n:
        return [Part(pos, ref, alt)]  # plain anchored insertion or deletion
    k = min(m, n)
    parts = [Part(pos + i, ref[i], alt[i]) for i in range(k) if ref[i] != alt[i]]
    if n > m:
        parts.append(Part(pos + m - 1, ref[m - 1], ref[m - 1] + alt[m:]))
    elif m > n:
        parts.append(Part(pos + n - 1, ref[n - 1 :], ref[n - 1]))
    return sorted(parts)


def apply_parts(sequence: str, parts: list[Part], offset: int = 1) -> str:
    """Apply parts (1-based, ``offset`` = genome position of ``sequence[0]``) to a sequence."""
    out: list[str] = []
    last = 0
    for part in sorted(parts, key=lambda p: (p.span0()[0], p.is_indel)):
        start, end = part.span0()
        start -= offset - 1
        end -= offset - 1
        if start < last:
            raise ValueError(f"overlapping parts at {part.pos}")
        out.append(sequence[last:start])
        if part.is_indel:
            out.append(part.alt[1:])
        else:
            out.append(part.alt)
        last = end
    out.append(sequence[last:])
    return "".join(out)


def conflicts(a: Part, b: Part) -> bool:
    """True if Mason cannot apply both parts on the same haplotype unambiguously."""
    a0, a1 = a.span0()
    b0, b1 = b.span0()
    a_point, b_point = a0 == a1, b0 == b1
    if a_point and b_point:
        return a0 == b0  # two insertions at the same point
    if a_point or b_point:
        p = a0 if a_point else b0
        s, e, other = (b0, b1, b) if a_point else (a0, a1, a)
        return s < p < e or (p == s and other.is_indel)
    return a0 < b1 and b0 < a1


@dataclass
class HaplotypeEvents:
    """Accepted parts of one haplotype: planted (any order) + background (position order)."""

    planted: list[Part]
    bg_max_end: int = -1
    bg_last_indel_start: int = -1

    def conflict_with(self, part: Part) -> Part | None:
        for other in self.planted:
            if conflicts(part, other):
                return other
        return None

    def background_conflict(self, part: Part) -> bool:
        """Background parts arrive sorted; compare with the planted parts and earlier ones."""
        if self.conflict_with(part) is not None:
            return True
        s, e = part.span0()
        if s == e:
            return s < self.bg_max_end or s == self.bg_last_indel_start
        if s < self.bg_max_end:
            return True
        return part.is_indel and s == self.bg_last_indel_start

    def accept_background(self, part: Part) -> None:
        s, e = part.span0()
        self.bg_max_end = max(self.bg_max_end, e)
        if part.is_indel:
            self.bg_last_indel_start = s
