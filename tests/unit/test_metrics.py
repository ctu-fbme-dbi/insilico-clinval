"""Metrics, acceptance criteria, regression diff and report (no external tools)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pysam
import pytest
import yaml
from pydantic import ValidationError
from typer.testing import CliRunner

from insilico_clinval.cli import app
from insilico_clinval.diff.run import concordance, diff
from insilico_clinval.fidelity.stats import wilson
from insilico_clinval.jsonio import write_json
from insilico_clinval.metrics.acceptance import evaluate
from insilico_clinval.metrics.compute import compute, criterion_code
from insilico_clinval.metrics.load import contig_group, happy_counts, load_run
from insilico_clinval.schemas.enums import Classification, EvaluationLevel, OutcomeClass
from insilico_clinval.schemas.fidelity import CallingMetric
from insilico_clinval.schemas.metrics import AcceptanceConfig, Criterion, MetricsReport
from insilico_clinval.schemas.outcome import (
    ExpectedOutcome,
    ExpectedOutcomeSet,
    OutcomeRecord,
    OutcomeSet,
)

P, LP, B, VUS = (
    Classification.PATHOGENIC,
    Classification.LIKELY_PATHOGENIC,
    Classification.BENIGN,
    Classification.UNCERTAIN_SIGNIFICANCE,
)
PAR = [("chrX", 10000, 2781479)]
C = OutcomeClass

# variant id -> (expected class, reportable, acceptable classes, outcome, observed class,
#                 reported, held_out)
PLANTED: dict[str, tuple[Any, ...]] = {
    "v_plp_ok": (P, True, (P, LP), C.CORRECT, P, True, "false"),
    "v_plp_vus": (LP, True, (P, LP), C.MISCLASSIFIED, VUS, False, "true"),
    "v_plp_lost": (P, True, (P, LP), C.NOT_CALLED, None, None, "true"),
    "v_carrier": (P, False, (P, LP), C.REPORT_MISMATCH, P, True, "false"),
    "v_blb_ok": (B, False, (B,), C.CORRECT, B, False, "unknown"),
    "v_vus_ok": (VUS, False, (VUS,), C.CORRECT, VUS, False, "false"),
}


def _record(vid: str, spec: tuple[Any, ...], sample: str, pos: int) -> OutcomeRecord:
    _, reportable, _, outcome, observed, reported, _ = spec
    return OutcomeRecord(
        variant_id=vid, sample_id=sample, pipeline="p", chrom="chr1", pos=pos, ref="A", alt="G",
        evaluation_level=EvaluationLevel.BOTH, outcome_class=outcome, stage=None,
        expected_reportable=reportable, observed_classification=observed,
        observed_reported=reported,
        evidence={"finding_found": observed is not None,
                  "acmg_criteria": "PVS1,PM2,PP3" if observed is not None else None},
    )  # fmt: skip


def _false_report(sample: str, pos: int, call: str, chrom: str = "chrX") -> OutcomeRecord:
    return OutcomeRecord(
        variant_id=None, sample_id=sample, pipeline="p", chrom=chrom, pos=pos, ref="C", alt="T",
        evaluation_level=None, outcome_class=C.FALSE_REPORT, stage=None, observed_reported=True,
        evidence={"call": call, "planted_negative": False},
    )  # fmt: skip


def _happy_vcf(path: Path, rows: list[str]) -> Path:
    header = [
        "##fileformat=VCFv4.2", "##contig=<ID=chr1,length=10000000>",
        "##contig=<ID=chrX,length=156040895>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">',
        '##FORMAT=<ID=BD,Number=1,Type=String,Description="x">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY",
    ]  # fmt: skip
    plain = path.with_suffix("")
    plain.write_text("\n".join(header + rows) + "\n")
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return path


def _run(
    root: Path,
    pipeline: str = "p",
    *,
    changes: dict[str, tuple[Any, ...]] | None = None,
    fp_rows: int = 1,
) -> Path:
    """A run directory with one compared sample of ``pipeline`` (``changes`` override PLANTED)."""
    run = root / "suiteA"
    planted = {**PLANTED, **(changes or {})}
    compare = run / "sut" / pipeline / "compare"
    compare.mkdir(parents=True, exist_ok=True)
    expected = [
        ExpectedOutcome(
            variant_id=v, expected_classification=s[0], acceptable_classifications=list(s[2]),
            expected_reportable=s[1], reportable_basis="t", condition=None,
            evaluation_level=EvaluationLevel.BOTH,
            expected_acmg_criteria=["PVS1", "PM2_Supporting"],
        )
        for v, s in planted.items()
    ]  # fmt: skip
    write_json(run / "s.expected_outcomes.json",
               ExpectedOutcomeSet(schema_version=1, sample_id="s", outcomes=expected))  # fmt: skip
    records = [_record(v, s, "s", 100 * (i + 1)) for i, (v, s) in enumerate(planted.items())]
    reports = [_false_report("s", 5_000_000, "fp_call")] + (
        [_false_report("s", 6_000_000, "fp_call")] if fp_rows > 1 else []
    )
    write_json(compare / "s.outcomes.json", OutcomeSet(
        schema_version=1, sample_id="s", pipeline=pipeline, parameters={},
        records=records, false_reports=reports))  # fmt: skip
    (run / "held_out.tsv").write_text(
        "variant_id\theld_out\tapproved_after_cutoff\n"
        + "".join(f"{v}\t{'true' if s[6] in ('true', 'unknown') else 'false'}\t"
                  f"{'false' if s[6] == 'unknown' else ''}\n" for v, s in planted.items())
    )  # fmt: skip
    (run / "variants.selected.json").write_text(json.dumps(
        {"variants": [{"variant_id": v, "variant_type": "snv"} for v in planted]}))  # fmt: skip
    m = CallingMetric(dataset="s", type="SNP", stratum="*", subtype="*", recall=wilson(9, 10),
                      precision=wilson(9, 9), f1=0.9)  # fmt: skip
    (compare / "s.calling.json").write_text(json.dumps([m.model_dump(mode="json")]))
    happy = compare / "s.happy"
    happy.mkdir(exist_ok=True)
    rows = ["chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP",
            "chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT:BD\t0/1:FN\t./.:.",
            "chrX\t20000\t.\tC\tT\t.\tPASS\t.\tGT:BD\t./.:.\t0/1:FP"]  # fmt: skip
    rows += [
        f"chrX\t{5_000_000 + i}\t.\tC\tT\t.\tPASS\t.\tGT:BD\t./.:.\t0/1:FP" for i in range(fp_rows)
    ]
    _happy_vcf(happy / "s.vcf.gz", rows)
    return run


def _value(report: MetricsReport, metric: str, **subset: str) -> tuple[int, int]:
    row = next(r for r in report.proportions if r.metric == metric and r.subset == subset)
    return row.value.k, row.value.n


# -- metrics ----------------------------------------------------------------------------


def test_metrics(tmp_path: Path) -> None:
    report = compute("p", [load_run(_run(tmp_path), "p")], PAR)
    assert _value(report, "report_sensitivity") == (1, 3)
    assert _value(report, "report_sensitivity", held_out="true") == (0, 2)
    assert _value(report, "report_specificity") == (2, 3)  # carrier reported
    assert _value(report, "report_specificity", group="recessive_carrier") == (0, 1)
    assert _value(report, "report_specificity", group="B/LB", held_out="unknown") == (1, 1)
    assert _value(report, "class_exact") == (4, 5)
    assert _value(report, "class_tier", group="P/LP") == (1, 2)
    assert _value(report, "criterion_recall") == (10, 10)  # PM2_Supporting == PM2
    assert _value(report, "criterion_precision") == (10, 15)
    assert _value(report, "criterion_recall", criterion="PP3") == (0, 0)
    assert _value(report, "planted_calling") == (5, 6)
    row = next(r for r in report.proportions if r.metric == "report_specificity" and not r.subset)
    assert row.value == wilson(2, 3)
    counts = {(c.metric, tuple(sorted(c.subset.items()))): c.count for c in report.counts}
    assert counts[("false_positive_calls", (("contig_group", "chrX_PAR"),))] == 1
    assert counts[("false_positive_calls", (("contig_group", "chrX_nonPAR"),))] == 1
    assert counts[("false_reports", (("call", "fp_call"), ("planted_negative", "false")))] == 1
    assert _value(report, "calling_recall", contig_group="autosomes") == (1, 2)
    assert criterion_code("ps4_supporting") == "PS4"


def test_calling_metrics_are_summed(tmp_path: Path) -> None:
    a = load_run(_run(tmp_path / "a"), "p")
    b = load_run(_run(tmp_path / "b"), "p")
    report = compute("p", [a, b], PAR)
    assert _value(report, "calling_recall", type="SNP", stratum="*") == (18, 20)
    assert _value(report, "report_sensitivity", suite="suiteA") == (2, 6)


def test_contig_groups(tmp_path: Path) -> None:
    assert contig_group("chr7", 1, PAR) == "autosomes"
    assert contig_group("chrX", 20000, PAR) == "chrX_PAR"
    assert contig_group("chrX", 3_000_000, PAR) == "chrX_nonPAR"
    assert contig_group("chrX", 3_000_000, []) == "chrX"
    assert contig_group("chrM", 1, PAR) == "other"
    run = _run(tmp_path, fp_rows=3)
    counts = happy_counts(run / "sut/p/compare/s.happy/s.vcf.gz", PAR)
    assert counts[("chrX_nonPAR", "query_FP")] == 3 and counts[("autosomes", "truth_FN")] == 1


# -- acceptance -----------------------------------------------------------------------------


def test_acceptance(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="unknown metric"):
        Criterion(metric="nope", min_point=0.5)
    with pytest.raises(ValidationError, match="takes max_count only"):
        Criterion(metric="false_reports", min_point=0.5)
    with pytest.raises(ValidationError, match="set min_point"):
        Criterion(metric="report_sensitivity")
    report = compute("p", [load_run(_run(tmp_path), "p")], PAR)
    cfg = AcceptanceConfig(name="t", criteria=[
        Criterion(metric="report_specificity", min_point=0.6),
        Criterion(metric="report_sensitivity", min_lower=0.5),
        Criterion(metric="false_reports", subset={"call": "fp_call"}, max_count=0),
        Criterion(metric="report_sensitivity", subset={"held_out": "nope"}, min_point=0.5),
    ])  # fmt: skip
    assert [r.status for r in evaluate(cfg, report)] == ["pass", "fail", "fail", "not_evaluable"]


def test_metrics_cli_and_fail_on_reject(tmp_path: Path) -> None:
    run = _run(tmp_path)
    pipe = tmp_path / "p.yaml"
    pipe.write_text(yaml.safe_dump({
        "schema_version": 1, "name": "p",
        "secondary": {"adapter": "precomputed",
                      "precomputed": {"version": "1", "samples": {"s": {"raw_vcf": "x.vcf.gz"}}}},
        "tertiary": {"adapter": "none"},
    }))  # fmt: skip
    acc = tmp_path / "acc.yaml"
    acc.write_text(yaml.safe_dump(
        {"name": "t", "criteria": [{"metric": "report_sensitivity", "min_point": 0.9}]}
    ))  # fmt: skip
    args = ["metrics", "--pipeline", str(pipe), "--run", str(run), "--out", str(tmp_path / "m"),
            "--acceptance", str(acc)]  # fmt: skip
    runner = CliRunner()
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    first = (tmp_path / "m/metrics.json").read_bytes()
    assert runner.invoke(app, [*args, "--fail-on-reject"]).exit_code == 7
    assert (tmp_path / "m/metrics.json").read_bytes() == first  # deterministic
    tsv = (tmp_path / "m/metrics.tsv").read_text().splitlines()
    assert tsv[0].startswith("metric\tsubset\tk\tn")
    manifest = json.loads((tmp_path / "m/run_manifest.json").read_text())
    assert {"metrics_json", "metrics_tsv"} <= {o["role"] for o in manifest["outputs"]}
    md = tmp_path / "r.md"
    assert (
        runner.invoke(app, ["report", "--metrics", str(tmp_path / "m"), "--out", str(md)]).exit_code
        == 0
    )
    assert "Sensitivity" in md.read_text() and "fail" in md.read_text()


# -- diff -------------------------------------------------------------------------------


def test_diff_kinds_and_noise(tmp_path: Path) -> None:
    root = tmp_path
    _run(root, "base")
    _run(root, "cand", fp_rows=2, changes={
        "v_plp_ok": (P, True, (P, LP), C.NOT_CALLED, None, None, "false"),          # new failure
        "v_plp_vus": (LP, True, (P, LP), C.CORRECT, LP, True, "true"),              # recovered
        "v_plp_lost": (P, True, (P, LP), C.FILTERED, None, None, "true"),           # stage change
        "v_blb_ok": (B, False, (B,), C.CORRECT, B, False, "unknown"),
        "v_vus_ok": (VUS, False, (VUS,), C.CORRECT, VUS, False, "false"),
    })  # fmt: skip
    _run(root, "noise", changes={
        "v_plp_lost": (P, True, (P, LP), C.FILTERED, None, None, "true"),
    })  # fmt: skip
    report = diff([root / "suiteA"], "base", "cand", noise=None, par=PAR)
    kinds = {(c.kind, c.variant_id) for c in report.changes}
    assert kinds == {
        ("new_failure", "v_plp_ok"), ("recovered", "v_plp_vus"), ("stage_change", "v_plp_lost"),
        ("false_report_added", None),
    }  # fmt: skip
    assert report.verdict == "regression"
    noisy = diff([root / "suiteA"], "base", "cand", noise="noise", par=PAR)
    assert {c.variant_id for c in noisy.changes if c.noise} == {"v_plp_lost"}
    fp = {d.contig_group: d.counts for d in report.calling}["chrX_nonPAR"]
    assert (fp["baseline_query_FP"], fp["candidate_query_FP"]) == (1, 2)
    assert report.concordance == []  # no secondary.json: nothing to compare


def test_concordance(tmp_path: Path) -> None:
    def vcf(name: str, rows: list[str]) -> Path:
        plain = tmp_path / name
        plain.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n"
                         "##contig=<ID=chrX,length=156040895>\n"
                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">\n'
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                         + "\n".join(rows) + "\n")  # fmt: skip
        pysam.tabix_index(str(plain), preset="vcf", force=True)
        return Path(f"{plain}.gz")

    a = vcf("a.vcf", ["chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
                      "chr1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
                      "chrX\t3000000\t.\tC\tT\t.\tPASS\t.\tGT\t0/1"])  # fmt: skip
    b = vcf("b.vcf", ["chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
                      "chr1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t1/1",
                      "chrX\t3000000\t.\tC\tT\t.\tRefCall\t.\tGT\t0/0"])  # fmt: skip
    rows = {r.contig_group: (r.identical, r.only_baseline, r.only_candidate)
            for r in concordance("run", "s", a, b, PAR)}  # fmt: skip
    assert rows == {"autosomes": (1, 1, 1), "chrX_nonPAR": (0, 1, 0)}
