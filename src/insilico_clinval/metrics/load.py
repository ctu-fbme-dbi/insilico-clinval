"""Reading the outputs of `clinval compare` for one pipeline in one run directory."""

from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from insilico_clinval.schemas.enums import Classification, EvaluationLevel
from insilico_clinval.schemas.fidelity import CallingMetric
from insilico_clinval.schemas.finding import SecondaryResult
from insilico_clinval.schemas.outcome import (
    ExpectedOutcome,
    ExpectedOutcomeSet,
    OutcomeRecord,
    OutcomeSet,
)
from insilico_clinval.variants.held_out import HELD_OUT_TSV
from insilico_clinval.variants.report import SELECTED_JSON

P_LP = {Classification.PATHOGENIC, Classification.LIKELY_PATHOGENIC}
B_LB = {Classification.BENIGN, Classification.LIKELY_BENIGN}
Interval = tuple[str, int, int]  # 0-based, half-open


class MetricsInputError(ValueError):
    """A comparison output needed for the metrics is missing."""


def group_of(expected: ExpectedOutcome) -> str:
    """P/LP (reportable), recessive carrier (P/LP, not reportable), B/LB, VUS or calling-only."""
    cls = expected.expected_classification
    if expected.evaluation_level == EvaluationLevel.CALLING or cls is None:
        return "calling_only"
    if cls in P_LP:
        return "P/LP" if expected.expected_reportable else "recessive_carrier"
    if cls in B_LB:
        return "B/LB"
    return "VUS"


def tier(cls: Classification | None) -> str | None:
    if cls is None:
        return None
    return "P/LP" if cls in P_LP else "B/LB" if cls in B_LB else "VUS"


@dataclass(frozen=True)
class PlantedRow:
    run: str
    record: OutcomeRecord
    expected: ExpectedOutcome
    group: str
    held_out: str  # true, false, unknown
    variant_type: str  # snv, indel


@dataclass
class RunData:
    run: str
    planted: list[PlantedRow] = field(default_factory=list)
    false_reports: list[OutcomeRecord] = field(default_factory=list)
    calling: list[CallingMetric] = field(default_factory=list)
    happy_vcfs: dict[str, Path] = field(default_factory=dict)
    raw_vcfs: dict[str, Path] = field(default_factory=dict)


def _held_out_labels(run_dir: Path) -> dict[str, str]:
    """held_out per variant; held out but approved before the snapshot = ``unknown``."""
    path = run_dir / HELD_OUT_TSV
    if not path.is_file():
        return {}
    labels: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ambiguous = r["held_out"] == "true" and r.get("approved_after_cutoff") == "false"
            labels[r["variant_id"]] = "unknown" if ambiguous else r["held_out"]
    return labels


def _variant_types(run_dir: Path) -> dict[str, str]:
    path = run_dir / SELECTED_JSON
    if not path.is_file():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {
        v["variant_id"]: "snv" if v["variant_type"] == "snv" else "indel" for v in data["variants"]
    }


def load_run(run_dir: Path, pipeline: str) -> RunData:
    """All compared samples of ``pipeline`` in ``run_dir`` (sorted by sample)."""
    compare = run_dir / "sut" / pipeline / "compare"
    sets = sorted(compare.glob("*.outcomes.json"))
    if not sets:
        raise MetricsInputError(f"{compare}: no outcomes (run clinval compare first)")
    held_out, types = _held_out_labels(run_dir), _variant_types(run_dir)
    data = RunData(run=str(run_dir))
    for path in sets:
        s = OutcomeSet.model_validate_json(path.read_text(encoding="utf-8"))
        exp_path = run_dir / f"{s.sample_id}.expected_outcomes.json"
        if not exp_path.is_file():
            raise MetricsInputError(f"{exp_path} not found")
        expected = {
            o.variant_id: o
            for o in ExpectedOutcomeSet.model_validate_json(exp_path.read_text()).outcomes
        }
        for r in s.records:
            assert r.variant_id is not None
            e = expected[r.variant_id]
            data.planted.append(
                PlantedRow(
                    run=data.run,
                    record=r,
                    expected=e,
                    group=group_of(e),
                    held_out=held_out.get(r.variant_id, "unknown"),
                    variant_type=types.get(
                        r.variant_id, "snv" if len(r.ref) == len(r.alt) else "indel"
                    ),
                )
            )
        data.false_reports.extend(s.false_reports)
        calling = compare / f"{s.sample_id}.calling.json"
        if calling.is_file():
            data.calling += [
                CallingMetric.model_validate(m) for m in json.loads(calling.read_text())
            ]
        happy = compare / f"{s.sample_id}.happy" / f"{s.sample_id}.vcf.gz"
        if happy.is_file():
            data.happy_vcfs[s.sample_id] = happy
        sec = run_dir / "sut" / pipeline / f"{s.sample_id}.secondary.json"
        if sec.is_file():
            raw = SecondaryResult.model_validate_json(sec.read_text()).raw_vcf
            data.raw_vcfs[s.sample_id] = Path(raw) if Path(raw).is_absolute() else run_dir / raw
    return data


# -- contig groups -----------------------------------------------------------------------------


def read_par(path: Path | None) -> list[Interval]:
    if path is None:
        return []
    out: list[Interval] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        if len(f) >= 3 and not line.startswith(("#", "track")):
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def contig_group(chrom: str, pos: int, par: Sequence[Interval]) -> str:
    """autosomes, chrX/chrY (split into _PAR and _nonPAR when PAR regions are given), other."""
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return "autosomes"
    if name not in ("X", "Y"):
        return "other"
    sex = f"chr{name}"
    if not par:
        return sex
    in_par = any(c == chrom and s < pos <= e for c, s, e in par)
    return f"{sex}_PAR" if in_par else f"{sex}_nonPAR"


def happy_counts(path: Path, par: Sequence[Interval]) -> Counter[tuple[str, str]]:
    """(contig group, truth_TP | truth_FN | query_FP) counts of a hap.py output VCF."""
    counts: Counter[tuple[str, str]] = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            group = contig_group(rec.chrom, rec.pos, par)
            truth, query = rec.samples["TRUTH"].get("BD"), rec.samples["QUERY"].get("BD")
            if truth in ("TP", "FN"):
                counts[(group, f"truth_{truth}")] += 1
            if query == "FP":
                counts[(group, "query_FP")] += 1
    return counts
