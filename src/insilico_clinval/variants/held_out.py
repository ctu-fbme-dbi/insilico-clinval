"""Held-out status of the eRepo variants of an existing run, recomputed from a reference.

The status is normally set at selection (``sources.clingen_erepo.held_out_reference``). Adding a
reference to a suite afterwards would change its configuration and invalidate the run directory,
so ``clinval held-out`` writes ``held_out.tsv`` next to the selection instead; the comparison picks
it up. As a plausibility check, each allele is also compared with the approval date of its eRepo
records: an allele whose records were all approved after the reference was taken should be held
out, and one approved before it should be in the reference. The two disagree for alleles that
were re-approved after the snapshot, that the reference failed to convert, or that changed
representation, so disagreement is reported, not treated as an error.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from insilico_clinval.genome.alleles import left_normalize
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.variant import TestVariant
from insilico_clinval.variants.erepo import HeldOutReference

HELD_OUT_TSV = "held_out.tsv"
COLUMNS = (
    "variant_id", "gene", "erepo_id", "held_out", "first_approval", "approved_after_cutoff",
    "agrees_with_date",
)  # fmt: skip


@dataclass(frozen=True)
class HeldOutRow:
    variant_id: str
    gene: str | None
    erepo_id: str
    held_out: bool
    first_approval: str | None
    approved_after_cutoff: bool | None

    @property
    def agrees_with_date(self) -> bool | None:
        if self.approved_after_cutoff is None:
            return None
        return self.held_out == self.approved_after_cutoff


def annotate(
    variants: list[TestVariant],
    reference_set: HeldOutReference,
    reference: Reference,
    *,
    approvals: dict[str, str] | None = None,
    cutoff: str | None = None,
) -> list[HeldOutRow]:
    """One row per eRepo variant, sorted by variant id."""
    rows: list[HeldOutRow] = []
    for v in sorted(variants, key=lambda t: t.variant_id):
        if v.clingen is None:
            continue
        ids = [v.clingen.erepo_id, *v.clingen.other_erepo_ids, v.clingen.caid]
        pos, ref, alt = left_normalize(v.chrom, v.pos, v.ref, v.alt, reference)
        held_out = not reference_set.contains((v.chrom, pos, ref, alt), ids)
        dates = sorted(approvals[i] for i in ids if i is not None and approvals and i in approvals)
        first = dates[0] if dates else None
        after = None if first is None or cutoff is None else first > cutoff
        rows.append(HeldOutRow(v.variant_id, v.gene, v.clingen.erepo_id, held_out, first, after))
    return rows


def _cell(value: object) -> str:
    if value is None:
        return ""
    return ("true" if value else "false") if isinstance(value, bool) else str(value)


def write_tsv(rows: list[HeldOutRow], path: Path) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for r in rows:
            writer.writerow([_cell(getattr(r, c)) for c in COLUMNS])
    return path


def read_held_out(path: Path) -> dict[str, bool]:
    """variant_id -> held_out from a `held_out.tsv` (empty if the file does not exist)."""
    if not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8") as fh:
        return {
            r["variant_id"]: r["held_out"] == "true" for r in csv.DictReader(fh, delimiter="\t")
        }
