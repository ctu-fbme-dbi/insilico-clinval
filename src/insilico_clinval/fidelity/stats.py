"""Parsers of the raw tool outputs of the fidelity experiment and the Wilson interval.

All parsers take the text written by the tool, so the assembly step can re-read raw outputs
produced earlier (or on another machine) without running the tools again.
"""

from __future__ import annotations

import csv
import gzip
import math
from collections.abc import Iterator
from pathlib import Path

from insilico_clinval.schemas.fidelity import CallingMetric, CallsetQc, Proportion, ReadQc

Z95 = 1.959963984540054


class StatsParseError(ValueError):
    """A raw tool output lacks a required field."""


def wilson(k: int, n: int, z: float = Z95) -> Proportion:
    """k of n with the Wilson score interval (rounded to 6 decimals)."""
    if n == 0:
        return Proportion(k=k, n=n, value=None, lower=None, upper=None)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return Proportion(
        k=k,
        n=n,
        value=round(p, 6),
        lower=round(max(0.0, centre - half), 6),
        upper=round(min(1.0, centre + half), 6),
    )


# -- BED -------------------------------------------------------------------------------------


def merged_bed(path: Path) -> Iterator[tuple[str, int, int]]:
    """Merged half-open intervals of a BED (plain or gzip), contigs in first-seen order."""
    raw: dict[str, list[tuple[int, int]]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    for chrom, intervals in raw.items():
        intervals.sort()
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                yield chrom, cur_s, cur_e
                cur_s, cur_e = s, e
        yield chrom, cur_s, cur_e


def write_samtools_regions(bed: Path, out: Path) -> int:
    """Write `samtools stats -t` regions (1-based, inclusive); return the region length in bp."""
    total = 0
    with out.open("w", encoding="utf-8") as fh:
        for chrom, s, e in merged_bed(bed):
            fh.write(f"{chrom}\t{s + 1}\t{e}\n")
            total += e - s
    return total


def region_length(regions: Path) -> int:
    """Length in bp of a region file written by :func:`write_samtools_regions`."""
    total = 0
    with regions.open(encoding="utf-8") as fh:
        for line in fh:
            _, s, e = line.rstrip("\n").split("\t")
            total += int(e) - int(s) + 1
    return total


# -- samtools stats --------------------------------------------------------------------------


def _sn(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith("SN\t"):
            f = line.split("\t")
            out[f[1].rstrip(":")] = f[2]
    return out


def _need(sn: dict[str, str], key: str, source: str) -> float:
    if key not in sn:
        raise StatsParseError(f"{source}: missing SN field {key!r}")
    return float(sn[key])


def coverage_metrics(hist: dict[int, int], region_bp: int) -> tuple[float, float, float | None]:
    """Mean depth, fraction >= 20x and fold-80 over ``region_bp`` bases (uncovered bases = 0)."""
    covered = sum(hist.values())
    if region_bp <= 0 or covered > region_bp:
        raise StatsParseError(
            f"coverage histogram ({covered} bp) exceeds the regions ({region_bp} bp)"
        )
    depths = {0: region_bp - covered, **hist} if covered < region_bp else dict(hist)
    mean = sum(d * n for d, n in depths.items()) / region_bp
    ge20 = sum(n for d, n in depths.items() if d >= 20) / region_bp
    target = 0.2 * region_bp
    cumulative = 0
    p20 = 0
    for d in sorted(depths):
        cumulative += depths[d]
        if cumulative >= target:
            p20 = d
            break
    fold80 = round(mean / p20, 4) if p20 > 0 else None
    return round(mean, 4), round(ge20, 6), fold80


def parse_samtools_stats(text: str, *, dataset: str, region_bp: int) -> ReadQc:
    """ReadQc from `samtools stats` output (COV: depth histogram; last bin `[1000<]` = 1000)."""
    src = f"samtools stats ({dataset})"
    sn = _sn(text)
    total = _need(sn, "raw total sequences", src)
    bases = _need(sn, "bases mapped (cigar)", src)
    indels = 0
    hist: dict[int, int] = {}
    for line in text.splitlines():
        if line.startswith("IC\t"):
            indels += sum(int(x) for x in line.split("\t")[2:6])
        elif line.startswith("COV\t"):
            f = line.split("\t")
            hist[int(f[2])] = hist.get(int(f[2]), 0) + int(f[3])
    mean, ge20, fold80 = coverage_metrics(hist, region_bp)

    def frac(key: str) -> float:
        return round(_need(sn, key, src) / total, 6) if total else 0.0

    return ReadQc(
        dataset=dataset,
        reads_total=int(total),
        reads_mapped=frac("reads mapped"),
        reads_properly_paired=frac("reads properly paired"),
        reads_mq0=frac("reads MQ0"),
        insert_size_mean=_need(sn, "insert size average", src),
        insert_size_sd=_need(sn, "insert size standard deviation", src),
        error_rate=_need(sn, "error rate", src),
        indel_rate=round(indels / bases, 9) if bases else 0.0,
        coverage_mean=mean,
        coverage_ge20=ge20,
        fold80=fold80,
    )


# -- bcftools stats --------------------------------------------------------------------------


def parse_bcftools_stats(text: str, *, dataset: str) -> CallsetQc:
    """CallsetQc from `bcftools stats -s -` output (first sample of the PSC section)."""
    src = f"bcftools stats ({dataset})"
    counts: dict[str, int] = {}
    ti_tv: float | None = None
    het_hom: float | None = None
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "SN":
            counts[f[2].rstrip(":")] = int(f[3])
        elif f[0] == "TSTV" and ti_tv is None:
            ts, tv = int(f[2]), int(f[3])
            ti_tv = round(ts / tv, 4) if tv else None
        elif f[0] == "PSC" and het_hom is None:
            hom_alt, het = int(f[4]), int(f[5])
            het_hom = round(het / hom_alt, 4) if hom_alt else None
    for key in ("number of SNPs", "number of indels"):
        if key not in counts:
            raise StatsParseError(f"{src}: missing SN field {key!r}")
    return CallsetQc(
        dataset=dataset,
        snvs=counts["number of SNPs"],
        indels=counts["number of indels"],
        ti_tv=ti_tv,
        het_hom=het_hom,
    )


# -- hap.py ----------------------------------------------------------------------------------

HAPPY_COLUMNS = (
    "Type", "Subtype", "Subset", "Filter", "Genotype", "METRIC.F1_Score",
    "TRUTH.TP", "TRUTH.FN", "QUERY.TP", "QUERY.FP",
)  # fmt: skip


def _int(value: str) -> int:
    return int(float(value)) if value not in ("", "nan", ".") else 0


def parse_happy_extended(text: str, *, dataset: str, pass_only: bool) -> list[CallingMetric]:
    """CallingMetric rows of a hap.py `*.extended.csv` (genotype `*`, filter PASS or ALL)."""
    reader = csv.DictReader(text.splitlines())
    missing = [c for c in HAPPY_COLUMNS if c not in (reader.fieldnames or [])]
    if missing:
        raise StatsParseError(f"hap.py extended.csv ({dataset}): missing columns {missing}")
    wanted = "PASS" if pass_only else "ALL"
    rows: list[CallingMetric] = []
    for r in reader:
        if r["Type"] not in ("SNP", "INDEL") or r["Filter"] != wanted or r["Genotype"] != "*":
            continue
        if r.get("QQ", "*") not in ("*", ""):
            continue
        t_tp, t_fn = _int(r["TRUTH.TP"]), _int(r["TRUTH.FN"])
        q_tp, q_fp = _int(r["QUERY.TP"]), _int(r["QUERY.FP"])
        f1 = r["METRIC.F1_Score"]
        rows.append(
            CallingMetric(
                dataset=dataset,
                type=r["Type"],  # type: ignore[arg-type]
                stratum=r["Subset"],
                subtype=r["Subtype"],
                recall=wilson(t_tp, t_tp + t_fn),
                precision=wilson(q_tp, q_tp + q_fp),
                f1=round(float(f1), 6) if f1 not in ("", "nan") else None,
            )
        )
    rows.sort(key=lambda m: (m.type, m.stratum, m.subtype))
    return rows
