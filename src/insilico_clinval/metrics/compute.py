"""Metrics over the comparison outputs of one pipeline, with 95 % Wilson intervals.

Report level (primary): sensitivity on planted variants that must be reported, specificity on
planted variants that must not (B/LB, VUS, recessive carriers). A variant lost before the report
counts as not reported. Classification level (secondary): 5-tier, acceptable and 3-tier
concordance and per-criterion agreement with the curation, over variants the tool classified.
Calling: planted variants called with the right genotype, hap.py recall and precision summed over
samples, and false-positive calls per contig group. Failures and false reports are counts.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path

from insilico_clinval.fidelity.stats import wilson
from insilico_clinval.metrics.load import Interval, PlantedRow, RunData, happy_counts, tier
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import OutcomeClass
from insilico_clinval.schemas.metrics import CountRow, MetricRow, MetricsReport

CALLING_FAILURES = {
    OutcomeClass.NOT_SIMULATED, OutcomeClass.NOT_COVERED, OutcomeClass.MISMAPPED,
    OutcomeClass.NOT_CALLED, OutcomeClass.FILTERED, OutcomeClass.GENOTYPE_MISMATCH,
}  # fmt: skip


def criterion_code(code: str) -> str:
    """ACMG code without its strength modifier (PM2_Supporting -> PM2)."""
    return code.strip().split("_")[0].upper()


def suite_of(row: PlantedRow) -> str:
    return Path(row.run).name


def _dims(row: PlantedRow, dim: str) -> str:
    if dim == "suite":
        return suite_of(row)
    return str(getattr(row, dim))


def _proportions(
    metric: str,
    rows: Sequence[PlantedRow],
    success: Callable[[PlantedRow], bool],
    subsets: Iterable[tuple[str, ...]],
) -> list[MetricRow]:
    out: list[MetricRow] = []
    for dims in subsets:
        groups: dict[tuple[str, ...], list[PlantedRow]] = {}
        for r in rows:
            groups.setdefault(tuple(_dims(r, d) for d in dims), []).append(r)
        for values, members in groups.items():
            k = sum(success(r) for r in members)
            subset = dict(zip(dims, values, strict=True))
            out.append(MetricRow(metric=metric, subset=subset, value=wilson(k, len(members))))
    return out


def _criteria(rows: Sequence[PlantedRow]) -> list[MetricRow]:
    found = [r for r in rows if r.record.evidence.get("finding_found") is True]
    exp_n: Counter[str] = Counter()
    obs_n: Counter[str] = Counter()
    both: Counter[str] = Counter()
    for r in found:
        expected = {criterion_code(c) for c in r.expected.expected_acmg_criteria}
        raw = r.record.evidence.get("acmg_criteria")
        observed = {criterion_code(c) for c in str(raw or "").split(",") if c.strip()}
        exp_n.update(expected)
        obs_n.update(observed)
        both.update(expected & observed)
    hits = sum(both.values())
    out = [
        MetricRow(metric="criterion_recall", subset={}, value=wilson(hits, sum(exp_n.values()))),
        MetricRow(metric="criterion_precision", subset={}, value=wilson(hits, sum(obs_n.values()))),
    ]
    for code in sorted(set(exp_n) | set(obs_n)):
        sub = {"criterion": code}
        out.append(
            MetricRow(metric="criterion_recall", subset=sub, value=wilson(both[code], exp_n[code]))
        )
        out.append(
            MetricRow(
                metric="criterion_precision", subset=sub, value=wilson(both[code], obs_n[code])
            )
        )
    return out


def _calling(runs: Sequence[RunData]) -> list[MetricRow]:
    sums: dict[tuple[str, str], list[int]] = {}
    for run in runs:
        for m in run.calling:
            if m.subtype != "*":
                continue
            s = sums.setdefault((m.type, m.stratum), [0, 0, 0, 0])
            s[0] += m.recall.k
            s[1] += m.recall.n
            s[2] += m.precision.k
            s[3] += m.precision.n
    out: list[MetricRow] = []
    for (vtype, stratum), (rk, rn, pk, pn) in sums.items():
        subset = {"type": vtype, "stratum": stratum}
        out.append(MetricRow(metric="calling_recall", subset=subset, value=wilson(rk, rn)))
        out.append(MetricRow(metric="calling_precision", subset=subset, value=wilson(pk, pn)))
    return out


def _sort_key(row: MetricRow | CountRow) -> tuple[str, tuple[tuple[str, str], ...]]:
    return row.metric, tuple(sorted(row.subset.items()))


def compute(pipeline: str, runs: Sequence[RunData], par: Sequence[Interval]) -> MetricsReport:
    planted = [r for run in runs for r in run.planted]
    interpreted = [r for r in planted if r.group != "calling_only"]
    reportable = [r for r in interpreted if r.expected.expected_reportable]
    negative = [r for r in interpreted if not r.expected.expected_reportable]
    classified = [r for r in interpreted if r.record.observed_classification is not None]

    def reported(r: PlantedRow) -> bool:
        return r.record.observed_reported is True

    def not_reported(r: PlantedRow) -> bool:
        return not reported(r)

    def exact(r: PlantedRow) -> bool:
        return r.record.observed_classification == r.expected.expected_classification

    def acceptable(r: PlantedRow) -> bool:
        return r.record.observed_classification in r.expected.acceptable_classifications

    def same_tier(r: PlantedRow) -> bool:
        observed, expected = r.record.observed_classification, r.expected.expected_classification
        return tier(observed) == tier(expected)

    def called(r: PlantedRow) -> bool:
        return r.record.outcome_class not in CALLING_FAILURES

    by_group: list[tuple[str, ...]] = [(), ("group",)]
    rows: list[MetricRow] = []
    rows += _proportions(
        "report_sensitivity", reportable, reported, [(), ("held_out",), ("suite",)]
    )
    rows += _proportions(
        "report_specificity",
        negative,
        not_reported,
        [(), ("group",), ("held_out",), ("group", "held_out"), ("suite",)],
    )
    rows += _proportions("class_exact", classified, exact, by_group)
    rows += _proportions("class_acceptable", classified, acceptable, by_group)
    rows += _proportions("class_tier", classified, same_tier, by_group)
    rows += _criteria(interpreted)
    rows += _proportions("planted_calling", planted, called, [(), ("variant_type",), ("suite",)])
    rows += _calling(runs)

    counts: list[CountRow] = []
    stage = Counter((r.record.outcome_class.value, r.group) for r in planted)
    for (cls, group), n in stage.items():
        counts.append(
            CountRow(metric="stage_counts", subset={"outcome_class": cls, "group": group}, count=n)
        )
    calls = Counter(
        (str(f.evidence.get("call")), str(f.evidence.get("planted_negative")).lower())
        for run in runs
        for f in run.false_reports
    )
    for (call, neg), n in calls.items():
        subset = {"call": call, "planted_negative": neg}
        counts.append(CountRow(metric="false_reports", subset=subset, count=n))
    per_sample = Counter(
        (Path(run.run).name, f.sample_id) for run in runs for f in run.false_reports
    )
    for (suite, sample), n in per_sample.items():
        subset = {"suite": suite, "sample_id": sample}
        counts.append(CountRow(metric="false_reports", subset=subset, count=n))

    contig: Counter[tuple[str, str]] = Counter()
    for run in runs:
        for vcf in run.happy_vcfs.values():
            contig.update(happy_counts(vcf, par))
    for group in sorted({g for g, _ in contig}):
        fp = contig[(group, "query_FP")]
        counts.append(
            CountRow(metric="false_positive_calls", subset={"contig_group": group}, count=fp)
        )
        tp, fn = contig[(group, "truth_TP")], contig[(group, "truth_FN")]
        rows.append(
            MetricRow(
                metric="calling_recall", subset={"contig_group": group}, value=wilson(tp, tp + fn)
            )
        )

    return MetricsReport(
        schema_version=SCHEMA_VERSION,
        pipeline=pipeline,
        runs=[run.run for run in runs],
        proportions=sorted(rows, key=_sort_key),
        counts=sorted(counts, key=_sort_key),
    )
