"""Planting completeness across run directories (Figure 4A)."""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

from insilico_clinval.fidelity.stats import wilson
from insilico_clinval.schemas.enums import PlantingStatus, Zygosity
from insilico_clinval.schemas.fidelity import PlantingRow, PlantingSummary, ZygositySummary
from insilico_clinval.schemas.simulation import PlantingCheck

PLANTING_COLUMNS = (
    "run", "sample_id", "variant_id", "zygosity", "expected_vaf", "spanning_reads", "alt_reads",
    "vaf", "status",
)  # fmt: skip


class PlantingError(ValueError):
    """A run directory has no planting check."""


def collect_planting(run_dirs: list[Path]) -> list[PlantingRow]:
    """One row per planted variant and sample, sorted by run, sample and variant."""
    rows: list[PlantingRow] = []
    for run_dir in run_dirs:
        checks = sorted(run_dir.glob("*.planting_check.json"))
        if not checks:
            raise PlantingError(f"{run_dir}: no *.planting_check.json (run `clinval simulate`)")
        for path in checks:
            check = PlantingCheck.model_validate_json(path.read_text(encoding="utf-8"))
            rows += [
                PlantingRow(
                    run=run_dir.name,
                    sample_id=check.sample_id,
                    variant_id=v.variant_id,
                    zygosity=v.zygosity,
                    expected_vaf=v.expected_vaf,
                    spanning_reads=v.spanning_reads,
                    alt_reads=v.alt_reads,
                    vaf=v.vaf,
                    status=v.status,
                )
                for v in check.variants
            ]
    rows.sort(key=lambda r: (r.run, r.sample_id, r.variant_id))
    return rows


def summarise_planting(rows: list[PlantingRow]) -> PlantingSummary:
    by_zygosity: list[ZygositySummary] = []
    for zyg in Zygosity:
        sel = [r for r in rows if r.zygosity == zyg]
        if not sel:
            continue
        vafs = [r.vaf for r in sel if r.vaf is not None]
        alts = [r.alt_reads for r in sel]
        by_zygosity.append(
            ZygositySummary(
                zygosity=zyg,
                variants=len(sel),
                vaf_median=round(statistics.median(vafs), 4) if vafs else None,
                vaf_min=min(vafs) if vafs else None,
                vaf_max=max(vafs) if vafs else None,
                alt_reads_median=float(statistics.median(alts)),
                alt_reads_min=min(alts),
            )
        )
    ok = sum(r.status == PlantingStatus.OK for r in rows)
    return PlantingSummary(
        runs=sorted({r.run for r in rows}),
        samples=len({(r.run, r.sample_id) for r in rows}),
        variants=len(rows),
        complete=wilson(ok, len(rows)),
        by_status={s.value: n for s in PlantingStatus if (n := sum(r.status == s for r in rows))},
        by_zygosity=by_zygosity,
    )


def write_planting_tsv(rows: list[PlantingRow], path: Path) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PLANTING_COLUMNS)
        for r in rows:
            w.writerow(
                [
                    r.run, r.sample_id, r.variant_id, r.zygosity.value, r.expected_vaf,
                    r.spanning_reads, r.alt_reads, "" if r.vaf is None else r.vaf, r.status.value,
                ]
            )  # fmt: skip
    return path


def read_planting_tsv(path: Path) -> list[PlantingRow]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [
            PlantingRow(
                run=r["run"],
                sample_id=r["sample_id"],
                variant_id=r["variant_id"],
                zygosity=Zygosity(r["zygosity"]),
                expected_vaf=float(r["expected_vaf"]),
                spanning_reads=int(r["spanning_reads"]),
                alt_reads=int(r["alt_reads"]),
                vaf=float(r["vaf"]) if r["vaf"] else None,
                status=PlantingStatus(r["status"]),
            )
            for r in csv.DictReader(fh, delimiter="\t")
        ]
