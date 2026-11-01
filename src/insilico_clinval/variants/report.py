"""Selection outputs: variants.selected.json, variants.selected.tsv, selection_report.md."""

from __future__ import annotations

import csv
import io
from pathlib import Path

from insilico_clinval.config import LoadedSuite
from insilico_clinval.jsonio import write_json
from insilico_clinval.schemas.enums import SampleRole, VariantCategory
from insilico_clinval.schemas.selection import SelectionResult
from insilico_clinval.variants.select import selection_failures, variants_per_sample

SELECTED_JSON = "variants.selected.json"
SELECTED_TSV = "variants.selected.tsv"
SELECTION_REPORT = "selection_report.md"

TSV_COLUMNS = [
    "variant_id",
    "rule_id",
    "category",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_type",
    "indel_length",
    "gene",
    "source_database",
    "source_release",
    "source_record_id",
    "rsid",
    "clnsig",
    "review_status",
    "stars",
    "vcep",
    "erepo_classification",
    "met_criteria",
    "held_out",
    "strata",
    "selection_reason",
]


def selected_tsv(result: SelectionResult) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_COLUMNS)
    for v in result.variants:
        writer.writerow(
            [
                v.variant_id,
                v.rule_id,
                v.category.value,
                v.chrom,
                v.pos,
                v.ref,
                v.alt,
                v.variant_type.value,
                v.indel_length,
                v.gene or "",
                v.source.database,
                v.source.release,
                v.source.record_id,
                v.source.rsid or "",
                "|".join(v.clinvar.clnsig) if v.clinvar else "",
                v.clinvar.review_status if v.clinvar else "",
                v.clinvar.stars if v.clinvar else "",
                v.clingen.vcep if v.clingen else "",
                v.clingen.classification.value if v.clingen else "",
                ",".join(v.clingen.met_criteria) if v.clingen else "",
                str(v.held_out).lower(),
                ",".join(v.strata),
                v.selection_reason,
            ]
        )
    return buf.getvalue()


# Roy et al. 2018 (AMP/CAP): at least 59 variants per category for 95 % reliability at 95 %
# confidence.
MIN_PER_CATEGORY = 59


def _pairing_section(suite: LoadedSuite, result: SelectionResult) -> list[str]:
    """Per-gene counts of curated variants in positive vs negative samples (with held-out)."""
    roles = {s.sample_id: s.role for s in suite.config.samples}
    per_sample = variants_per_sample(result)
    counts: dict[str, dict[str, int]] = {}
    vceps: dict[str, set[str]] = {}
    for sample_id, variants in per_sample.items():
        role = roles[sample_id]
        if role == SampleRole.MIXED:
            continue
        for v in variants:
            if v.category != VariantCategory.CURATED or not v.gene:
                continue
            c = counts.setdefault(v.gene, {})
            c[role.value] = c.get(role.value, 0) + 1
            if v.held_out:
                c[f"{role.value}_held_out"] = c.get(f"{role.value}_held_out", 0) + 1
            if v.clingen is not None:
                vceps.setdefault(v.gene, set()).add(v.clingen.vcep)
    lines = ["## Positive / negative pairing (curated variants per gene)", ""]
    if not counts:
        return [*lines, "No curated variants in positive or negative samples.", ""]
    lines += [
        "| gene | positive (P/LP) | of which held out | negative (B/LB) | of which held out "
        "| VCEP | note |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    totals = {"positive": 0, "negative": 0}
    for gene in sorted(counts):
        c = counts[gene]
        pos, neg = c.get("positive", 0), c.get("negative", 0)
        totals["positive"] += pos
        totals["negative"] += neg
        note = "" if pos == neg else "**asymmetric**"
        vcep = "; ".join(sorted(vceps.get(gene, set())))
        lines.append(
            f"| {gene} | {pos} | {c.get('positive_held_out', 0)} | {neg} | "
            f"{c.get('negative_held_out', 0)} | {vcep} | {note} |"
        )
    lines.append("")
    for group, n in totals.items():
        if n < MIN_PER_CATEGORY:
            lines.append(
                f"- {group}: {n} curated variants, below the {MIN_PER_CATEGORY} per category "
                "recommended by Roy et al. 2018 (95 % reliability, 95 % confidence)"
            )
    return [*lines, ""]


def selection_report(suite: LoadedSuite, result: SelectionResult) -> str:
    cfg = suite.config
    lines = [
        f"# Selection report — {cfg.suite.name} {cfg.suite.version}",
        "",
        f"- seed: `{cfg.suite.seed}`",
        f"- reference: {cfg.reference.build} (`{Path(cfg.reference.fasta).name}`)",
    ]
    erepo = cfg.sources.clingen_erepo
    if erepo is not None:
        ref = erepo.held_out_reference
        held = f"held-out reference {ref.release}" if ref else "no held-out reference"
        lines.append(f"- ClinGen eRepo release: {erepo.release} ({held})")
    if cfg.sources.clinvar is not None:
        lines.append(f"- ClinVar release: {cfg.sources.clinvar.release}")
    if cfg.sources.giab is not None:
        lines.append(f"- {cfg.sources.giab.name} release: {cfg.sources.giab.release}")
    if cfg.background is not None:
        bg = cfg.background
        with_bg = [s.sample_id for s in cfg.samples if cfg.uses_background(s)]
        lines.append(
            f"- background: {bg.name} {bg.release} ({bg.sex.value}) in samples "
            f"{', '.join(with_bg) or 'none'}"
        )
    lines += [
        f"- min_distance_bp: {cfg.constraints.min_distance_bp}; "
        f"max_indel_length: {cfg.constraints.max_indel_length}",
        f"- selected variants: {len(result.variants)}",
        "",
    ]
    failures = selection_failures(result, [s.sample_id for s in cfg.samples])
    if failures:
        lines += ["## FAILED", ""] + [f"- {f}" for f in failures] + [""]
    else:
        lines += ["## Status: OK", ""]

    for rule in result.rules:
        status = "satisfied" if rule.satisfied else "UNSATISFIED"
        req = "required" if rule.required else "optional"
        lines += [
            f"## Rule `{rule.rule_id}` — {status} ({req})",
            "",
            f"Samples: {', '.join(rule.samples)}",
            "",
            "| filter | candidates remaining |",
            "|---|---:|",
        ]
        lines += [f"| {name} | {n} |" for name, n in rule.filter_counts]
        lines += ["", "| bucket | requested | selected |", "|---|---:|---:|"]
        lines += [f"| {b} | {n} | {rule.selected.get(b, 0)} |" for b, n in rule.requested.items()]
        lines.append("")
        if rule.missing:
            lines += ["Missing:", ""] + [f"- {m}" for m in rule.missing] + [""]
        if rule.rejections:
            lines += ["Rejected by deduplication:", ""]
            lines += [f"- `{r.variant_id}`: {r.reason}" for r in rule.rejections]
            lines.append("")

    lines += _pairing_section(suite, result)
    per_sample = variants_per_sample(result)
    lines += ["## Samples", "", "| sample | role | sex | variants |", "|---|---|---|---:|"]
    for s in cfg.samples:
        n = len(per_sample.get(s.sample_id, []))
        lines.append(f"| {s.sample_id} | {s.role.value} | {s.sex.value} | {n} |")
    lines.append("")
    if result.warnings:
        lines += ["## Warnings", ""] + [f"- {w}" for w in result.warnings] + [""]
    return "\n".join(lines)


def write_selection_outputs(
    suite: LoadedSuite, result: SelectionResult, out_dir: Path
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "variants_selected_json": write_json(out_dir / SELECTED_JSON, result),
        "variants_selected_tsv": out_dir / SELECTED_TSV,
        "selection_report": out_dir / SELECTION_REPORT,
    }
    paths["variants_selected_tsv"].write_text(selected_tsv(result), encoding="utf-8")
    paths["selection_report"].write_text(selection_report(suite, result), encoding="utf-8")
    return paths
