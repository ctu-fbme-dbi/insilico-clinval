"""Compare the outputs of one pipeline with the truth, sample by sample.

Layout (``<run>/sut/<pipeline>/compare/``)::

    <sample>.regions.bed          scored regions (see compare.regions)
    <sample>.planted.bed          planted-variant windows (stratum ``planted``)
    <sample>.happy/strata.tsv     hap.py stratification (planted windows + optional strata)
    <sample>.happy/<sample>.*     raw hap.py output
    <sample>.calling.json         CallingMetric rows (PASS calls)
    <sample>.outcomes.json        OutcomeSet
    outcomes.tsv                  all OutcomeRecords of all compared samples
"""

from __future__ import annotations

import csv
import json
import logging
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pysam

from insilico_clinval._runner import Toolbox
from insilico_clinval.adapters import LoadedPipeline
from insilico_clinval.adapters.run import load_secondary, sut_dir
from insilico_clinval.compare.attribute import VariantEvidence, attribute, evidence_dict
from insilico_clinval.compare.evidence import (
    Allele,
    FindingIndex,
    alignment_counts,
    calls_near,
    happy_calls,
    normalized,
    open_alignment,
    pick,
)
from insilico_clinval.compare.regions import (
    Interval,
    evaluation_regions,
    planted_window,
    read_segments,
    split,
    write_bed,
)
from insilico_clinval.fidelity.stats import parse_happy_extended
from insilico_clinval.genome.reference import Reference
from insilico_clinval.jsonio import write_json
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import EvaluationLevel, OutcomeClass, Sex
from insilico_clinval.schemas.finding import CanonicalFinding, FindingSet
from insilico_clinval.schemas.outcome import ExpectedOutcomeSet, OutcomeRecord, OutcomeSet
from insilico_clinval.schemas.simulation import PlantingCheck
from insilico_clinval.simulate.check import variant_window
from insilico_clinval.variants.held_out import HELD_OUT_TSV, read_held_out

log = logging.getLogger(__name__)


class CompareError(ValueError):
    """A required input of the comparison is missing or inconsistent."""


@dataclass(frozen=True)
class CompareParams:
    window_pad: int = 50
    mismapped_fraction: float = 0.5
    edge_bp: int = 500
    threads: int = 4
    happy_gender: bool = False  # pass the sample's sex to hap.py (--gender)

    def as_dict(self) -> dict[str, float | int | str]:
        return {
            "window_pad": self.window_pad,
            "mismapped_fraction": self.mismapped_fraction,
            "edge_bp": self.edge_bp,
            "happy_gender": "sample" if self.happy_gender else "auto",
        }


@dataclass(frozen=True)
class HappyJob:
    truth: Path
    query: Path
    reference: Path
    regions: Path
    strata: Path
    prefix: Path
    threads: int
    mounts: tuple[Path, ...] = ()
    gender: Sex | None = None  # None: hap.py guesses (treats chrX as diploid)


HappyRunner = Callable[[HappyJob, Toolbox, ManifestBuilder | None, str], None]


def happy_args(job: HappyJob) -> list[str]:
    return [
        str(job.truth), str(job.query), "-r", str(job.reference), "-f", str(job.regions),
        "-o", str(job.prefix), "--engine", "vcfeval", "--engine-vcfeval-path", "rtg",
        "--threads", str(job.threads), "--pass-only", "--stratification", str(job.strata),
        *(["--gender", job.gender.value] if job.gender is not None else []),
    ]  # fmt: skip


def run_happy(job: HappyJob, toolbox: Toolbox, manifest: ManifestBuilder | None, step: str) -> None:
    job.prefix.parent.mkdir(parents=True, exist_ok=True)
    mounts = [job.truth.parent, job.query.parent, job.reference.parent, job.regions.parent]
    mounts += [job.prefix.parent, *job.mounts]
    toolbox.run(
        "hap.py",
        happy_args(job),
        step=step,
        manifest=manifest,
        mounts=mounts,
        cwd=job.prefix.parent,
        stderr_path=job.prefix.parent / "happy.stderr.log",
    )


def compare_dir(run_dir: Path, pipeline: LoadedPipeline) -> Path:
    return sut_dir(run_dir, pipeline) / "compare"


def _load(path: Path, what: str) -> str:
    if not path.is_file():
        raise CompareError(f"{what} not found: {path}")
    return path.read_text(encoding="utf-8")


def _rel_to(run_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    p = Path(value)
    return p if p.is_absolute() else run_dir / p


def write_strata(out: Path, strata: Sequence[tuple[str, Path]]) -> Path:
    """hap.py stratification TSV (name, absolute BED path)."""
    path = out / "strata.tsv"
    path.write_text("".join(f"{n}\t{p.resolve()}\n" for n, p in strata), encoding="utf-8")
    return path


def compare_sample(
    run_dir: Path,
    pipeline: LoadedPipeline,
    sample_id: str,
    *,
    reference: Reference,
    strata: Sequence[tuple[str, Path]],
    toolbox: Toolbox,
    manifest: ManifestBuilder | None,
    params: CompareParams,
    happy: HappyRunner = run_happy,
    sex: Sex | None = None,
) -> OutcomeSet:
    out = compare_dir(run_dir, pipeline)
    out.mkdir(parents=True, exist_ok=True)
    name = pipeline.config.name
    expected = ExpectedOutcomeSet.model_validate_json(
        _load(run_dir / f"{sample_id}.expected_outcomes.json", "expected outcomes")
    )
    planting = PlantingCheck.model_validate_json(
        _load(run_dir / f"{sample_id}.planting_check.json", "planting check")
    )
    simulation = json.loads(_load(run_dir / f"{sample_id}.simulation.json", "simulation result"))
    if "segments" not in simulation:
        raise CompareError(f"{sample_id}.simulation.json names no segments file")
    secondary = load_secondary(run_dir, pipeline, sample_id)
    truth = run_dir / f"{sample_id}.truth.vcf.gz"
    raw_vcf = _rel_to(run_dir, secondary.raw_vcf)
    query = _rel_to(run_dir, secondary.filtered_vcf) or raw_vcf
    alignment = _rel_to(run_dir, secondary.alignment)
    assert raw_vcf is not None and query is not None
    findings: list[CanonicalFinding] = []
    tertiary = pipeline.config.tertiary.adapter != "none"
    if tertiary:
        findings = FindingSet.model_validate_json(
            _load(sut_dir(run_dir, pipeline) / f"{sample_id}.findings.json", "findings")
        ).findings

    planted = {v.variant_id: v for v in planting.variants}
    missing = [o.variant_id for o in expected.outcomes if o.variant_id not in planted]
    if missing:
        raise CompareError(f"{sample_id}: planted variants without planting check: {missing}")

    # regions and hap.py
    windows = [
        planted_window(v.chrom, v.pos, v.ref, v.alt, params.window_pad) for v in planting.variants
    ]
    planted_bed = write_bed(sorted(windows, key=lambda w: reference.sort_key(w[0], w[1])),
                            out / f"{sample_id}.planted.bed")  # fmt: skip
    scored = evaluation_regions(
        read_segments(run_dir / simulation["segments"]), windows, reference, edge_bp=params.edge_bp
    )
    regions = write_bed(split(scored), out / f"{sample_id}.regions.bed")
    happy_dir = out / f"{sample_id}.happy"
    happy_dir.mkdir(parents=True, exist_ok=True)
    strata_tsv = write_strata(happy_dir, [("planted", planted_bed.resolve()), *strata])
    prefix = happy_dir / sample_id
    job = HappyJob(
        truth=truth.resolve(),
        query=query.resolve(),
        reference=reference.fasta.resolve(),
        regions=regions.resolve(),
        strata=strata_tsv.resolve(),
        prefix=prefix.resolve(),
        threads=params.threads,
        mounts=tuple(p.resolve().parent for _, p in strata),
        gender=sex if params.happy_gender else None,
    )
    happy(job, toolbox, manifest, f"compare hap.py ({name}, {sample_id})")
    calling = parse_happy_extended(
        _load(Path(f"{prefix}.extended.csv"), "hap.py extended.csv"),
        dataset=sample_id,
        pass_only=True,
    )
    write_json(out / f"{sample_id}.calling.json", calling)

    # per planted variant
    held_out = read_held_out(run_dir / HELD_OUT_TSV)
    index = FindingIndex(findings, reference)
    records: list[OutcomeRecord] = []
    planted_alleles: dict[Allele, tuple[str, bool, EvaluationLevel]] = {}
    split_windows: list[tuple[str, int, int, str]] = []  # found by hap.py, allele not annotated
    aln = open_alignment(alignment, reference.fasta) if alignment is not None else None
    with pysam.VariantFile(f"{prefix}.vcf.gz") as hvcf, pysam.VariantFile(str(raw_vcf)) as rvcf:
        for o in sorted(expected.outcomes, key=lambda e: e.variant_id):
            p = planted[o.variant_id]
            allele = normalized(reference, (p.chrom, p.pos, p.ref, p.alt))
            planted_alleles[allele] = (o.variant_id, o.expected_reportable, o.evaluation_level)
            s0, e0 = variant_window(reference, p.chrom, p.pos, p.ref, p.alt)
            counts = (
                alignment_counts(aln, reference, (p.chrom, p.pos, p.ref, p.alt)) if aln else None
            )
            hcall, how = pick(happy_calls(hvcf, reference, "TRUTH", p.chrom, s0, e0), allele)
            near = calls_near(rvcf, reference, p.chrom, s0, e0)
            nonpass = [c.filters for c in near if c.allele == allele and not c.passed]
            finding = index.find(allele) if tertiary else None
            ev = VariantEvidence(
                planting_status=p.status,
                true_alt_reads=p.alt_reads,
                min_alt_reads=planting.min_alt_reads,
                aligned_alt_reads=counts.alt if counts else None,
                happy_decision=hcall.decision if hcall else None,
                happy_kind=hcall.kind if hcall else None,
                nonpass_calls=nonpass,
                tertiary=tertiary,
                finding_found=finding is not None,
                observed_classification=finding.classification if finding else None,
                observed_reported=finding.reported if finding else None,
            )
            result = attribute(o, ev, mismapped_fraction=params.mismapped_fraction)
            if tertiary and finding is None and hcall is not None and hcall.decision == "TP":
                split_windows.append((p.chrom, s0, e0, o.variant_id))
            records.append(
                OutcomeRecord(
                    variant_id=o.variant_id,
                    sample_id=sample_id,
                    pipeline=name,
                    chrom=p.chrom,
                    pos=p.pos,
                    ref=p.ref,
                    alt=p.alt,
                    evaluation_level=o.evaluation_level,
                    outcome_class=result.outcome,
                    stage=result.stage,
                    expected_classification=o.expected_classification,
                    expected_reportable=o.expected_reportable,
                    observed_classification=ev.observed_classification,
                    observed_reported=ev.observed_reported,
                    gene=finding.gene if finding else None,
                    evidence=evidence_dict(
                        ev,
                        happy_match=how,
                        aligned_spanning=counts.spanning if counts else None,
                        aligned_alt_mapq0=counts.alt_mapq0 if counts else None,
                        calls_in_window=len(near),
                        acmg_criteria=",".join(finding.acmg_criteria) if finding else None,
                        held_out=held_out.get(o.variant_id),
                    ),
                )
            )
        false_reports, linked = _false_reports(
            findings,
            planted_alleles,
            hvcf,
            reference,
            sample_id=sample_id,
            pipeline=name,
            scored=scored,
            split_windows=split_windows,
        )
        for r in records:
            if r.variant_id in linked:
                r.evidence["split_report"] = ",".join(linked[r.variant_id])
    if aln is not None:
        aln.close()

    result_set = OutcomeSet(
        schema_version=SCHEMA_VERSION,
        sample_id=sample_id,
        pipeline=name,
        parameters=params.as_dict(),
        records=records,
        false_reports=false_reports,
    )
    path = write_json(out / f"{sample_id}.outcomes.json", result_set)
    if manifest is not None:
        manifest.add_output(f"outcomes:{name}:{sample_id}", path)
        manifest.add_output(f"calling:{name}:{sample_id}", out / f"{sample_id}.calling.json")
    log.info(
        "%s %s: %d planted, %d correct, %d false reports",
        name, sample_id, len(records),
        sum(r.outcome_class == OutcomeClass.CORRECT for r in records), len(false_reports),
    )  # fmt: skip
    return result_set


def _false_reports(
    findings: Sequence[CanonicalFinding],
    planted: dict[Allele, tuple[str, bool, EvaluationLevel]],
    hvcf: pysam.VariantFile,
    reference: Reference,
    *,
    sample_id: str,
    pipeline: str,
    scored: Sequence[Interval] = (),
    split_windows: Sequence[tuple[str, int, int, str]] = (),
) -> tuple[list[OutcomeRecord], dict[str, list[str]]]:
    """Reported findings that match no planted variant expected to be reported.

    A planted variant evaluated at ``calling`` level only has no expected report decision, so
    reporting it is neither correct nor a false report (its record is scored on calling).

    A reported finding inside the window of a planted variant that hap.py found but whose allele
    was not annotated is the same failure in another representation (e.g. a delins called as an
    insertion and an SNV): it is linked to that variant (returned per variant id) instead of being
    counted again as a false report.

    ``call`` of a false report: ``true_variant`` / ``fp_call`` from the exact hap.py match;
    otherwise ``outside_scored_regions`` (not in the scored segments and windows) or
    ``no_exact_match`` (scored, but hap.py has no record with exactly this allele).
    """
    out: list[OutcomeRecord] = []
    linked: dict[str, list[str]] = {}
    for f in findings:
        if not f.reported:
            continue
        allele = normalized(reference, (f.chrom, f.pos, f.ref, f.alt))
        variant_id, reportable, level = planted.get(allele, (None, False, None))
        if variant_id is not None and (reportable or level == EvaluationLevel.CALLING):
            continue
        s0, e0 = variant_window(reference, f.chrom, f.pos, f.ref, f.alt)
        if variant_id is None:
            split_of = next(
                (vid for c, a, b, vid in split_windows if c == f.chrom and s0 < b and a < e0),
                None,
            )
            if split_of is not None:
                linked.setdefault(split_of, []).append(f"{f.chrom}:{f.pos}:{f.ref}:{f.alt}")
                continue
        qcall, how = pick(happy_calls(hvcf, reference, "QUERY", f.chrom, s0, e0), allele)
        decision = qcall.decision if qcall is not None and how == "exact" else None
        call = {"TP": "true_variant", "FP": "fp_call"}.get(decision or "")
        if call is None:
            inside = any(c == f.chrom and a < e0 and s0 < b for c, a, b in scored)
            call = "no_exact_match" if inside else "outside_scored_regions"
        out.append(
            OutcomeRecord(
                variant_id=None,
                sample_id=sample_id,
                pipeline=pipeline,
                chrom=f.chrom,
                pos=f.pos,
                ref=f.ref,
                alt=f.alt,
                evaluation_level=None,
                outcome_class=OutcomeClass.FALSE_REPORT,
                stage=None,
                observed_classification=f.classification,
                observed_reported=True,
                gene=f.gene,
                evidence={
                    "call": call,
                    "happy_decision": decision,
                    "planted_negative": variant_id is not None,
                    "planted_variant_id": variant_id,
                    "acmg_criteria": ",".join(f.acmg_criteria),
                },
            )
        )
    out.sort(key=lambda r: (reference.sort_key(r.chrom, r.pos), r.ref, r.alt))
    return out, linked


TSV_COLUMNS = (
    "pipeline", "sample_id", "variant_id", "chrom", "pos", "ref", "alt", "evaluation_level",
    "outcome_class", "stage", "expected_classification", "expected_reportable",
    "observed_classification", "observed_reported", "gene", "call", "planted_negative",
    "held_out",
)  # fmt: skip


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(getattr(value, "value", value))


def write_outcomes_tsv(out: Path) -> Path:
    """`outcomes.tsv` from every `<sample>.outcomes.json` in ``out`` (sorted by sample)."""
    path = out / "outcomes.tsv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for f in sorted(out.glob("*.outcomes.json")):
            s = OutcomeSet.model_validate_json(f.read_text(encoding="utf-8"))
            for r in [*s.records, *s.false_reports]:
                row = r.model_dump()
                row["call"] = r.evidence.get("call")
                row["planted_negative"] = r.evidence.get("planted_negative")
                row["held_out"] = r.evidence.get("held_out")
                writer.writerow([_cell(row.get(c, getattr(r, c, None))) for c in TSV_COLUMNS])
    return path
