"""Acceptance criteria: loading and evaluation against a MetricsReport."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import ValidationError

from insilico_clinval.schemas.metrics import (
    COUNT_METRICS,
    AcceptanceConfig,
    Criterion,
    CriterionResult,
    MetricsReport,
)


class AcceptanceError(ValueError):
    """The acceptance configuration is missing or invalid."""


def load_acceptance(path: Path) -> AcceptanceConfig:
    if not path.is_file():
        raise AcceptanceError(f"acceptance config not found: {path}")
    try:
        return AcceptanceConfig.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    except (ValidationError, yaml.YAMLError) as exc:
        raise AcceptanceError(f"{path}: {exc}") from exc


def _judge(c: Criterion, report: MetricsReport) -> CriterionResult:
    if c.metric in COUNT_METRICS:
        # counts are split by several keys: sum every row whose subset contains the criterion's
        rows = [
            r
            for r in report.counts
            if r.metric == c.metric and c.subset.items() <= r.subset.items()
        ]
        total = sum(r.count for r in rows)
        assert c.max_count is not None
        passed = total <= c.max_count
        return CriterionResult(
            criterion=c, status="pass" if passed else "fail", observed=str(total)
        )
    match = [r for r in report.proportions if r.metric == c.metric and r.subset == c.subset]
    if not match or match[0].value.n == 0:
        return CriterionResult(criterion=c, status="not_evaluable", observed="no data")
    v = match[0].value
    assert v.value is not None and v.lower is not None
    ok = (c.min_point is None or v.value >= c.min_point) and (
        c.min_lower is None or v.lower >= c.min_lower
    )
    observed = f"{v.k}/{v.n} = {v.value:.4f} (95% CI {v.lower:.4f}-{v.upper:.4f})"
    return CriterionResult(criterion=c, status="pass" if ok else "fail", observed=observed)


def evaluate(config: AcceptanceConfig, report: MetricsReport) -> list[CriterionResult]:
    return [_judge(c, report) for c in config.criteria]
