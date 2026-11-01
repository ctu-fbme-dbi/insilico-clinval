"""Regression diff between a baseline and a candidate pipeline on the same simulated samples.

Planted variants are matched by (run, sample, variant id) and each change gets one kind:
``new_failure`` (correct -> failed), ``recovered`` (failed -> correct), ``stage_change`` (failed
in both, at another stage) or ``class_change`` (same correctness, but the observed class or the
report decision differs). False reports are matched by (run, sample, allele) and are ``added`` or
``removed``. With a noise run (a second run of the baseline), every change also seen between the
baseline and the noise run is labelled ``noise`` and does not count towards the verdict.

Calling is compared per contig group from the hap.py output VCFs; when the two pipelines have
different secondary outputs, their PASS calls are also compared record by record (run-to-run
concordance).
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

import pysam

from insilico_clinval.metrics.load import Interval, RunData, contig_group, happy_counts, load_run
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import OutcomeClass
from insilico_clinval.schemas.metrics import CallingDelta, Change, Concordance, DiffReport
from insilico_clinval.schemas.outcome import OutcomeRecord

Key = tuple[str, str, str]  # run, sample, variant id or allele


def _state(r: OutcomeRecord) -> str:
    cls = r.observed_classification.value if r.observed_classification else "-"
    rep = "-" if r.observed_reported is None else str(r.observed_reported).lower()
    return f"{r.outcome_class.value}/{cls}/{rep}"


def _allele(r: OutcomeRecord) -> str:
    return f"{r.chrom}:{r.pos}:{r.ref}>{r.alt}"


def _change(kind: str, run: str, r: OutcomeRecord, base: str, cand: str) -> Change:
    return Change.model_validate(
        {
            "kind": kind, "run": run, "sample_id": r.sample_id, "variant_id": r.variant_id,
            "chrom": r.chrom, "pos": r.pos, "ref": r.ref, "alt": r.alt, "gene": r.gene,
            "baseline": base, "candidate": cand, "noise": False,
        }
    )  # fmt: skip


def outcome_changes(a: RunData, b: RunData) -> list[Change]:
    """Changes from ``a`` (baseline) to ``b`` for one run directory."""
    out: list[Change] = []
    pa = {(r.record.sample_id, r.record.variant_id): r.record for r in a.planted}
    pb = {(r.record.sample_id, r.record.variant_id): r.record for r in b.planted}
    correct = OutcomeClass.CORRECT
    for key in sorted(set(pa) & set(pb), key=str):
        x, y = pa[key], pb[key]
        kind = None
        if x.outcome_class == correct and y.outcome_class != correct:
            kind = "new_failure"
        elif x.outcome_class != correct and y.outcome_class == correct:
            kind = "recovered"
        elif x.outcome_class != y.outcome_class:
            kind = "stage_change"
        elif (x.observed_classification, x.observed_reported) != (
            y.observed_classification,
            y.observed_reported,
        ):
            kind = "class_change"
        if kind is not None:
            out.append(_change(kind, a.run, y, _state(x), _state(y)))
    fa = {(f.sample_id, _allele(f)): f for f in a.false_reports}
    fb = {(f.sample_id, _allele(f)): f for f in b.false_reports}
    for key in sorted(set(fa) - set(fb)):
        f = fa[key]
        out.append(_change("false_report_removed", a.run, f, str(f.evidence.get("call")), "-"))
    for key in sorted(set(fb) - set(fa)):
        f = fb[key]
        out.append(_change("false_report_added", a.run, f, "-", str(f.evidence.get("call"))))
    return out


def _key(c: Change) -> Key:
    return c.run, c.sample_id, c.variant_id or f"{c.chrom}:{c.pos}:{c.ref}>{c.alt}"


def _pass_calls(path: Path, par: Sequence[Interval]) -> dict[str, set[tuple[object, ...]]]:
    """PASS records (contig, pos, ref, alts, genotype) of a VCF, by contig group."""
    out: dict[str, set[tuple[object, ...]]] = {}
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            if "PASS" not in rec.filter:
                continue
            gt = rec.samples[0].get("GT") if len(rec.samples) else None
            key = (rec.chrom, rec.pos, rec.ref, rec.alts, tuple(gt or ()))
            out.setdefault(contig_group(rec.chrom, rec.pos, par), set()).add(key)
    return out


def concordance(
    run: str, sample: str, a: Path, b: Path, par: Sequence[Interval]
) -> list[Concordance]:
    ca, cb = _pass_calls(a, par), _pass_calls(b, par)
    rows: list[Concordance] = []
    for group in sorted(set(ca) | set(cb)):
        x, y = ca.get(group, set()), cb.get(group, set())
        rows.append(
            Concordance(run=run, sample_id=sample, contig_group=group, identical=len(x & y),
                        only_baseline=len(x - y), only_candidate=len(y - x))
        )  # fmt: skip
    return rows


def calling_deltas(a: RunData, b: RunData, par: Sequence[Interval]) -> list[CallingDelta]:
    rows: list[CallingDelta] = []
    for sample in sorted(set(a.happy_vcfs) & set(b.happy_vcfs)):
        ha, hb = happy_counts(a.happy_vcfs[sample], par), happy_counts(b.happy_vcfs[sample], par)
        for group in sorted({g for g, _ in ha} | {g for g, _ in hb}):
            counts = {
                f"{who}_{what}": h[(group, what)]
                for who, h in (("baseline", ha), ("candidate", hb))
                for what in ("truth_TP", "truth_FN", "query_FP")
            }
            rows.append(
                CallingDelta(run=a.run, sample_id=sample, contig_group=group, counts=counts)
            )
    return rows


def diff(
    runs: Sequence[Path],
    baseline: str,
    candidate: str,
    *,
    noise: str | None,
    par: Sequence[Interval],
) -> DiffReport:
    changes: list[Change] = []
    calling: list[CallingDelta] = []
    conc: list[Concordance] = []
    for run_dir in runs:
        a, b = load_run(run_dir, baseline), load_run(run_dir, candidate)
        run_changes = outcome_changes(a, b)
        if noise is not None:
            noisy = {_key(c) for c in outcome_changes(a, load_run(run_dir, noise))}
            run_changes = [c.model_copy(update={"noise": _key(c) in noisy}) for c in run_changes]
        changes += run_changes
        calling += calling_deltas(a, b, par)
        for sample in sorted(set(a.raw_vcfs) & set(b.raw_vcfs)):
            ra, rb = a.raw_vcfs[sample], b.raw_vcfs[sample]
            if ra.resolve() != rb.resolve():
                conc += concordance(a.run, sample, ra, rb, par)
    counted = [c for c in changes if not c.noise]
    regression = any(c.kind in ("new_failure", "false_report_added") for c in counted)
    return DiffReport(
        schema_version=SCHEMA_VERSION,
        baseline=baseline,
        candidate=candidate,
        noise=noise,
        runs=[str(r) for r in runs],
        changes=changes,
        calling=calling,
        concordance=conc,
        verdict="regression" if regression else "no_regression",
    )


DIFF_COLUMNS = (
    "kind", "noise", "run", "sample_id", "variant_id", "chrom", "pos", "ref", "alt", "gene",
    "baseline", "candidate",
)  # fmt: skip


def _cell(value: object) -> str:
    if value is None:
        return ""
    return ("true" if value else "false") if isinstance(value, bool) else str(value)


def write_diff_tsv(report: DiffReport, path: Path) -> Path:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(DIFF_COLUMNS)
        for c in report.changes:
            row = c.model_dump()
            writer.writerow([_cell(row[k]) for k in DIFF_COLUMNS])
    return path
