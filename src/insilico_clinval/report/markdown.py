"""Markdown summaries of a MetricsReport and a DiffReport, and the metrics TSV."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from insilico_clinval.schemas.fidelity import Proportion
from insilico_clinval.schemas.metrics import DiffReport, MetricRow, MetricsReport


def row(*cells: object) -> str:
    """One Markdown table row."""
    return "| " + " | ".join(str(c) for c in cells) + " |"


def fmt(p: Proportion) -> str:
    """k/n = value % (lower to upper), percentages with one decimal."""
    if p.value is None:
        return f"{p.k}/{p.n} (n/a)"
    lo, hi = 100 * (p.lower or 0), 100 * (p.upper or 0)
    return f"{p.k}/{p.n} = {100 * p.value:.1f} % ({lo:.1f} to {hi:.1f})"


def _subset(s: dict[str, str]) -> str:
    return ", ".join(f"{k}={v}" for k, v in sorted(s.items())) or "all"


def write_metrics_tsv(report: MetricsReport, path: Path) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "subset", "k", "n", "value", "lower", "upper", "count"])
        for r in report.proportions:
            v = r.value
            w.writerow([r.metric, _subset(r.subset), v.k, v.n, v.value, v.lower, v.upper, ""])
        for c in report.counts:
            w.writerow([c.metric, _subset(c.subset), "", "", "", "", "", c.count])
    return path


def _table(rows: list[MetricRow], title: str) -> list[str]:
    if not rows:
        return []
    out = [f"### {title}", "", "| Subset | Value (95 % Wilson CI) |", "|---|---|"]
    out += [f"| {_subset(r.subset)} | {fmt(r.value)} |" for r in rows]
    return [*out, ""]


def metrics_markdown(report: MetricsReport) -> str:
    def pick(metric: str, keys: set[str] | None = None) -> list[MetricRow]:
        return [
            r for r in report.proportions
            if r.metric == metric and (keys is None or set(r.subset) <= keys)
        ]  # fmt: skip

    lines = [f"# Metrics — {report.pipeline}", "", "Runs: " + ", ".join(report.runs), ""]
    lines += ["## Report level (primary)", ""]
    lines += _table(pick("report_sensitivity"), "Sensitivity (P/LP reported)")
    lines += _table(pick("report_specificity"), "Specificity (B/LB, VUS, carriers not reported)")
    lines += ["## Classification (secondary)", ""]
    for m, title in (("class_exact", "5-tier exact"), ("class_acceptable", "Acceptable class"),
                     ("class_tier", "3-tier")):  # fmt: skip
        lines += _table(pick(m, {"group"}), title)
    lines += _table(pick("criterion_recall", set()), "ACMG criteria: recall (all codes)")
    lines += _table(pick("criterion_precision", set()), "ACMG criteria: precision (all codes)")
    lines += ["## Calling", ""]
    lines += _table(pick("planted_calling"), "Planted variants called with the right genotype")
    lines += _table(pick("calling_recall", {"contig_group"}), "Recall by contig group")
    lines += _table(pick("calling_recall", {"type", "stratum"}), "Recall by type and stratum")
    lines += _table(pick("calling_precision", {"type", "stratum"}), "Precision by type and stratum")
    lines += ["## Failures and false reports", "", "| Metric | Subset | Count |", "|---|---|---|"]
    lines += [f"| {c.metric} | {_subset(c.subset)} | {c.count} |" for c in report.counts]
    lines.append("")
    if report.criteria:
        lines += [f"## Acceptance ({report.acceptance})", "", "| Criterion | Status | Observed |",
                  "|---|---|---|"]  # fmt: skip
        for c in report.criteria:
            crit = c.criterion
            bounds = ", ".join(
                f"{k}={v}" for k, v in (("min_point", crit.min_point),
                                        ("min_lower", crit.min_lower),
                                        ("max_count", crit.max_count)) if v is not None
            )  # fmt: skip
            lines.append(
                f"| {crit.metric} [{_subset(crit.subset)}] {bounds} | {c.status} | {c.observed} |"
            )
        lines.append("")
    return "\n".join(lines)


def diff_markdown(report: DiffReport) -> str:
    kinds = Counter((c.kind, c.noise) for c in report.changes)
    lines = [
        f"# Regression diff — {report.candidate} vs {report.baseline}", "",
        f"Verdict: **{report.verdict}**"
        + (f" (noise run: {report.noise})" if report.noise else ""),
        "", "| Change | Counted | Noise |", "|---|---|---|",
    ]  # fmt: skip
    for kind in sorted({k for k, _ in kinds}):
        lines.append(f"| {kind} | {kinds[(kind, False)]} | {kinds[(kind, True)]} |")
    lines += ["", "## Changed items", ""]
    lines += [row("Kind", "Sample", "Variant / allele", "Gene", "Baseline", "Candidate")]
    lines += [row(*["---"] * 6)]
    for c in report.changes:
        name = c.variant_id or f"{c.chrom}:{c.pos}:{c.ref}>{c.alt}"
        label = c.kind + (" (noise)" if c.noise else "")
        lines.append(row(label, c.sample_id, name, c.gene or "", c.baseline, c.candidate))
    lines += ["", "## Calling by contig group", ""]
    lines += [row("Sample", "Group", "FP baseline", "FP candidate", "FN baseline", "FN candidate")]
    lines += [row(*["---"] * 6)]
    for d in report.calling:
        k = d.counts
        lines.append(
            row(d.sample_id, d.contig_group, k["baseline_query_FP"], k["candidate_query_FP"],
                k["baseline_truth_FN"], k["candidate_truth_FN"])
        )  # fmt: skip
    if report.concordance:
        lines += ["", "## Secondary concordance (PASS calls)", ""]
        lines += [row("Sample", "Group", "Identical", "Only baseline", "Only candidate")]
        lines += [row(*["---"] * 5)]
        for conc in report.concordance:
            lines.append(
                row(conc.sample_id, conc.contig_group, conc.identical, conc.only_baseline,
                    conc.only_candidate)
            )  # fmt: skip
    return "\n".join(lines) + "\n"
