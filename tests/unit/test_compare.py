"""Comparison with the truth and failure-stage attribution (no external tools)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pysam
import pytest

from insilico_clinval._runner import Toolbox
from insilico_clinval.adapters import load_pipeline
from insilico_clinval.compare.attribute import VariantEvidence, attribute
from insilico_clinval.compare.evidence import (
    FindingIndex,
    alignment_counts,
    calls_near,
    happy_calls,
    normalized,
    pick,
)
from insilico_clinval.compare.regions import evaluation_regions, planted_window, split
from insilico_clinval.compare.run import (
    CompareParams,
    HappyJob,
    compare_dir,
    compare_sample,
    happy_args,
    write_outcomes_tsv,
)
from insilico_clinval.genome.reference import Reference
from insilico_clinval.jsonio import write_json
from insilico_clinval.provenance import ManifestBuilder
from insilico_clinval.schemas.enums import (
    Classification,
    EvaluationLevel,
    OutcomeClass,
    PlantingStatus,
    Sex,
    Zygosity,
)
from insilico_clinval.schemas.finding import CanonicalFinding, FindingSet, SecondaryResult
from insilico_clinval.schemas.outcome import ExpectedOutcome, ExpectedOutcomeSet
from insilico_clinval.schemas.simulation import PlantingCheck, VariantPlanting
from tests.conftest import FIXTURES

P, LP, B = Classification.PATHOGENIC, Classification.LIKELY_PATHOGENIC, Classification.BENIGN
VUS = Classification.UNCERTAIN_SIGNIFICANCE


@pytest.fixture
def reference() -> Reference:
    return Reference(FIXTURES / "ref.fa")


def _base(reference: Reference, pos: int, chrom: str = "chr1") -> str:
    return reference.fetch(chrom, pos - 1, pos)


def _expected(
    variant_id: str = "v1",
    *,
    cls: Classification | None = P,
    acceptable: tuple[Classification, ...] = (P, LP),
    reportable: bool = True,
    level: EvaluationLevel = EvaluationLevel.BOTH,
) -> ExpectedOutcome:
    return ExpectedOutcome(
        variant_id=variant_id,
        expected_classification=cls,
        acceptable_classifications=list(acceptable),
        expected_reportable=reportable,
        reportable_basis="test",
        condition=None,
        evaluation_level=level,
    )


# -- regions ------------------------------------------------------------------------------


def test_evaluation_regions(reference: Reference) -> None:
    segments = [("chr1", 0, 10000), ("chrX", 2000, 8000)]  # chr2 not simulated
    windows = [
        planted_window("chrX", 2100, "A", "G", 50),  # in the trimmed segment end: still scored
        planted_window("chr1", 6000, "A", "G", 50),  # inside a segment: merged into it
    ]
    regions = evaluation_regions(segments, windows, reference, edge_bp=500)
    assert regions == [("chr1", 500, 9500), ("chrX", 2048, 2151), ("chrX", 2500, 7500)]
    assert split(regions, 4000) == [
        ("chr1", 500, 4500), ("chr1", 4500, 8500), ("chr1", 8500, 9500),
        ("chrX", 2048, 2151), ("chrX", 2500, 6500), ("chrX", 6500, 7500),
    ]  # fmt: skip


# -- attribution ----------------------------------------------------------------------

GOOD: dict[str, Any] = {
    "planting_status": PlantingStatus.OK,
    "true_alt_reads": 16,
    "min_alt_reads": 3,
    "aligned_alt_reads": 15,
    "happy_decision": "TP",
    "happy_kind": "gm",
    "finding_found": True,
    "observed_classification": P,
    "observed_reported": True,
}

NO_READS, NOT_ANNOTATED = PlantingStatus.NO_READS, OutcomeClass.NOT_ANNOTATED
CASES = [
    ({"planting_status": NO_READS, "true_alt_reads": 0}, OutcomeClass.NOT_SIMULATED, None),
    ({"true_alt_reads": 2, "aligned_alt_reads": 2}, OutcomeClass.NOT_COVERED, "simulated"),
    ({"aligned_alt_reads": 7, "happy_decision": "FN"}, OutcomeClass.MISMAPPED, "covered"),
    # a called allele is not mismapped, even if clipped reads hide it (long indels)
    ({"aligned_alt_reads": 2}, OutcomeClass.CORRECT, "reported"),
    ({"aligned_alt_reads": 2, "happy_decision": "FN", "happy_kind": "am"},
     OutcomeClass.GENOTYPE_MISMATCH, "passed"),
    ({"happy_decision": "FN", "happy_kind": "."}, OutcomeClass.NOT_CALLED, "mapped"),
    ({"happy_decision": "FN", "nonpass_calls": ["RefCall"]}, OutcomeClass.FILTERED, "called"),
    ({"happy_decision": "FN", "happy_kind": "am"}, OutcomeClass.GENOTYPE_MISMATCH, "passed"),
    ({"finding_found": False, "observed_classification": None}, NOT_ANNOTATED, "genotyped"),
    ({"observed_classification": VUS}, OutcomeClass.MISCLASSIFIED, "annotated"),
    ({"observed_reported": False}, OutcomeClass.REPORT_MISMATCH, "classified"),
    ({}, OutcomeClass.CORRECT, "reported"),
    ({"aligned_alt_reads": None}, OutcomeClass.CORRECT, "reported"),
]  # fmt: skip


@pytest.mark.parametrize(("change", "outcome", "stage"), CASES)
def test_attribution_classes(change: dict[str, Any], outcome: OutcomeClass, stage: str) -> None:
    ev = VariantEvidence(**{**GOOD, **change})
    result = attribute(_expected(), ev, mismapped_fraction=0.5)
    assert (result.outcome, result.stage) == (outcome, stage)


def test_attribution_calling_level_and_no_tertiary() -> None:
    ev = VariantEvidence(**{**GOOD, "finding_found": False, "observed_classification": None})
    calling = _expected(cls=None, acceptable=(), level=EvaluationLevel.CALLING)
    assert attribute(calling, ev, mismapped_fraction=0.5).outcome == OutcomeClass.CORRECT
    no_tertiary = VariantEvidence(**{**GOOD, "tertiary": False, "finding_found": False})
    result = attribute(_expected(), no_tertiary, mismapped_fraction=0.5)
    assert (result.outcome, result.stage) == (OutcomeClass.CORRECT, "genotyped")
    missed = VariantEvidence(**{**GOOD, "happy_decision": "FN", "happy_kind": "."})
    assert attribute(calling, missed, mismapped_fraction=0.5).outcome == OutcomeClass.NOT_CALLED


def test_attribution_negative_control() -> None:
    blb = _expected(cls=B, acceptable=(B, Classification.LIKELY_BENIGN), reportable=False)
    ok = VariantEvidence(**{**GOOD, "observed_classification": B, "observed_reported": False})
    assert attribute(blb, ok, mismapped_fraction=0.5).outcome == OutcomeClass.CORRECT
    wrong = VariantEvidence(**{**GOOD, "observed_classification": LP, "observed_reported": True})
    assert attribute(blb, wrong, mismapped_fraction=0.5).outcome == OutcomeClass.MISCLASSIFIED


# -- evidence ----------------------------------------------------------------------------------


def _write_vcf(path: Path, samples: list[str], rows: list[str]) -> Path:
    """bgzipped, indexed VCF with the given sample columns and data rows (FORMAT included)."""
    header = ["##fileformat=VCFv4.2", '##FILTER=<ID=RefCall,Description="x">']
    header += [f"##contig=<ID={c},length={n}>" for c, n in (("chr1", 10000), ("chrX", 8000))]
    header += [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="x">',
        '##FORMAT=<ID=BD,Number=1,Type=String,Description="x">',
        '##FORMAT=<ID=BK,Number=1,Type=String,Description="x">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]),
    ]  # fmt: skip
    header[-1] += "".join(f"\t{name}" for name in samples)
    plain = path.with_suffix("")
    plain.write_text("\n".join(header + rows) + "\n")
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return path


def test_calls_and_happy_lookup(tmp_path: Path, reference: Reference) -> None:
    b = _base(reference, 2001)
    alt = "A" if b != "A" else "C"
    d = reference.fetch("chr1", 2999, 3002)  # deletion of 2 bases after the anchor at 3000
    raw = _write_vcf(tmp_path / "raw.vcf.gz", ["S"], [
        f"chr1\t2001\t.\t{b}\t{alt}\t50\tPASS\t.\tGT\t0/1",
        f"chr1\t3000\t.\t{d}\t{d[0]}\t5\tRefCall\t.\tGT\t0/0",
    ])  # fmt: skip
    with pysam.VariantFile(str(raw)) as vcf:
        near = calls_near(vcf, reference, "chr1", 2990, 3010)
    assert [(c.allele, c.passed, c.filters) for c in near] == [
        (("chr1", 3000, d, d[0]), False, "RefCall")
    ]
    happy = _write_vcf(tmp_path / "happy.vcf.gz", ["TRUTH", "QUERY"], [
        f"chr1\t2001\t.\t{b}\t{alt}\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm",
        f"chr1\t3000\t.\t{d}\t{d[0]}\t.\tPASS\t.\tGT:BD:BK\t0/1:FN:am\t./.:.:.",
    ])  # fmt: skip
    with pysam.VariantFile(str(happy)) as vcf:
        truth = happy_calls(vcf, reference, "TRUTH", "chr1", 2990, 3010)
        query = happy_calls(vcf, reference, "QUERY", "chr1", 1990, 2010)
    call, how = pick(truth, ("chr1", 3000, d, d[0]))
    assert call is not None and (call.decision, call.kind, how) == ("FN", "am", "exact")
    assert [(c.decision, c.allele[1]) for c in query] == [("TP", 2001)]
    assert pick([], ("chr1", 1, "A", "C")) == (None, "none")


def test_finding_index_matches_other_representation(reference: Reference) -> None:
    ref3 = reference.fetch("chr1", 4999, 5003)  # 4 bases from 5000
    planted = normalized(reference, ("chr1", 5000, ref3[:3], ref3[0]))
    other = CanonicalFinding(
        sample_id="s", chrom="chr1", pos=5000, ref=ref3, alt=ref3[0] + ref3[3],
        classification=P, reported=True, tool="t", tool_version="v",
    )  # fmt: skip
    index = FindingIndex([other], reference)
    assert index.find(planted) is other
    assert index.find(("chr1", 7000, "A", "C")) is None


def test_alignment_counts(tmp_path: Path, reference: Reference) -> None:
    seq = reference.fetch("chr1", 0, 10000)
    pos = 2001
    alt = "A" if seq[pos - 1] != "A" else "C"
    path = tmp_path / "aln.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, start in enumerate(range(1940, 1990, 5)):
            read = seq[start : start + 100]
            if i % 2 == 0:
                read = read[: pos - 1 - start] + alt + read[pos - start :]
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = f"r{i}", 0, 0, start
            a.mapping_quality = 0 if i == 0 else 60
            a.cigarstring, a.query_sequence = "100M", read
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(a)
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path)) as aln:
        counts = alignment_counts(aln, reference, ("chr1", pos, seq[pos - 1], alt))
    assert (counts.spanning, counts.alt, counts.alt_mapq0) == (10, 5, 1)


# -- end to end with a recorded hap.py output -----------------------------------------


EXTENDED = (
    "Type,Subtype,Subset,Filter,Genotype,QQ,METRIC.F1_Score,TRUTH.TP,TRUTH.FN,QUERY.TP,QUERY.FP\n"
    "SNP,*,*,PASS,*,*,0.9,3,0,3,1\nINDEL,*,*,PASS,*,*,0.5,0,1,0,0\n"
    "SNP,*,planted,PASS,*,*,1.0,2,0,2,0\n"
)


def _run_dir(tmp_path: Path, reference: Reference) -> tuple[Path, Any, dict[str, Any]]:
    run = tmp_path / "run"
    run.mkdir()
    b = {p: _base(reference, p) for p in (1001, 2001, 3001, 4001, 7001, 8001, 9001, 9301)}
    dx = "A" if b[9001] != "A" else "C"  # delins: one base -> two, first base changed
    alt = {p: "A" if b[p] != "A" else "C" for p in b}
    d = reference.fetch("chr1", 5999, 6002)
    xb = _base(reference, 501, "chrX")
    xalt = "A" if xb != "A" else "C"
    variants = [
        ("v1_plp", "chr1", 1001, b[1001], alt[1001]),
        ("v2_plp_del", "chr1", 6000, d, d[0]),
        ("v3_blb", "chr1", 2001, b[2001], alt[2001]),
        ("v4_calling", "chr1", 8001, b[8001], alt[8001]),
        ("v5_delins", "chr1", 9001, b[9001], dx + "T"),
    ]
    planting = PlantingCheck(
        schema_version=1, sample_id="s", min_alt_reads=3, window_pad=3, complete=True,
        variants=[
            VariantPlanting(
                variant_id=v, chrom=c, pos=p, ref=r, alt=a, zygosity=Zygosity.HET,
                carrying_haplotypes=[1], expected_vaf=0.5, spanning_reads=30, alt_reads=15,
                ref_reads=15, other_reads=0, vaf=0.5, by_haplotype=[], status=PlantingStatus.OK,
            )
            for v, c, p, r, a in variants
        ],
    )  # fmt: skip
    write_json(run / "s.planting_check.json", planting)
    write_json(run / "s.expected_outcomes.json", ExpectedOutcomeSet(
        schema_version=1, sample_id="s",
        outcomes=[_expected("v1_plp"), _expected("v2_plp_del"),
                  _expected("v3_blb", cls=B, acceptable=(B,), reportable=False),
                  _expected("v4_calling", reportable=False, level=EvaluationLevel.CALLING),
                  _expected("v5_delins")],
    ))  # fmt: skip
    (run / "s.sim").mkdir()
    (run / "s.sim" / "segments.tsv").write_text(
        "segment\tchrom\tstart\tend\tploidy\tgroup\nchr1\tchr1\t1\t10000\t2\tdiploid\n"
    )
    (run / "s.simulation.json").write_text(json.dumps({"segments": "s.sim/segments.tsv"}))
    (run / "s.truth.vcf.gz").write_text("")
    (run / "held_out.tsv").write_text("variant_id\theld_out\nv1_plp\ttrue\nv3_blb\tfalse\n")
    raw = _write_vcf(run / "raw.vcf.gz", ["S"], [
        f"chr1\t1001\t.\t{b[1001]}\t{alt[1001]}\t50\tPASS\t.\tGT\t0/1",
        f"chr1\t2001\t.\t{b[2001]}\t{alt[2001]}\t50\tPASS\t.\tGT\t0/1",
        f"chr1\t4001\t.\t{b[4001]}\t{alt[4001]}\t50\tPASS\t.\tGT\t0/1",
        f"chr1\t7001\t.\t{b[7001]}\t{alt[7001]}\t50\tPASS\t.\tGT\t0/1",
        f"chr1\t8001\t.\t{b[8001]}\t{alt[8001]}\t50\tPASS\t.\tGT\t0/1",
    ])  # fmt: skip
    parser = {
        "format": "tsv",
        "columns": {"chrom": "c", "pos": "p", "ref": "r", "alt": "a", "classification": "k"},
        "classification_map": {"P": "pathogenic"},
    }
    command = {"tool_name": "fake", "version": "1", "runner": "local", "argv": ["true"]}
    cfg = {
        "schema_version": 1,
        "name": "p1",
        "secondary": {
            "adapter": "precomputed",
            "precomputed": {"version": "1", "samples": {"s": {"raw_vcf": "x.vcf.gz"}}},
        },
        "tertiary": {
            "adapter": "command",
            "command": {"command": command, "output": "out.tsv", "parser": parser},
        },
    }
    (tmp_path / "p.yaml").write_text(json.dumps(cfg))
    pipeline = load_pipeline(tmp_path / "p.yaml")
    sut = run / "sut" / "p1"
    sut.mkdir(parents=True)
    write_json(sut / "s.secondary.json", SecondaryResult(
        schema_version=1, sample_id="s", pipeline="p1", adapter="precomputed", tool="t",
        tool_version="1", alignment=None, raw_vcf=str(raw.relative_to(run)), filtered_vcf=None,
    ))  # fmt: skip

    def finding(
        pos: int, ref: str, alt_: str, cls: Classification, reported: bool, chrom: str = "chr1"
    ) -> CanonicalFinding:
        return CanonicalFinding(
            sample_id="s", chrom=chrom, pos=pos, ref=ref, alt=alt_, gene="G",
            classification=cls, reported=reported, tool="t", tool_version="1",
        )  # fmt: skip

    write_json(sut / "s.findings.json", FindingSet(
        schema_version=1, sample_id="s", pipeline="p1", adapter="command", tool="t",
        tool_version="1", input_vcf="x", findings=[
            finding(1001, b[1001], alt[1001], P, True),       # v1: correct
            finding(2001, b[2001], alt[2001], LP, True),      # v3: B/LB reported -> misclassified
            finding(4001, b[4001], alt[4001], P, True),       # background true variant
            finding(7001, b[7001], alt[7001], LP, True),      # false-positive call
            finding(8001, b[8001], alt[8001], P, True),       # v4 calling level: not a report
            finding(9001, b[9001], b[9001] + "T", P, True),   # v5 split: insertion part only
            finding(9301, b[9301], alt[9301], P, True),       # scored, no exact hap.py record
            finding(501, xb, xalt, P, True, chrom="chrX"),    # outside the scored regions
            finding(3001, b[3001], alt[3001], VUS, False),    # not reported: ignored
        ],
    ))  # fmt: skip
    happy_rows = [
        f"chr1\t1001\t.\t{b[1001]}\t{alt[1001]}\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm",
        f"chr1\t2001\t.\t{b[2001]}\t{alt[2001]}\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm",
        f"chr1\t4001\t.\t{b[4001]}\t{alt[4001]}\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm",
        f"chr1\t6000\t.\t{d}\t{d[0]}\t.\tPASS\t.\tGT:BD:BK\t0/1:FN:lm\t./.:.:.",
        f"chr1\t7001\t.\t{b[7001]}\t{alt[7001]}\t.\tPASS\t.\tGT:BD:BK\t./.:.:.\t0/1:FP:.",
        f"chr1\t8001\t.\t{b[8001]}\t{alt[8001]}\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm",
        f"chr1\t9001\t.\t{b[9001]}\t{dx}T\t.\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t./.:.:.",
    ]  # fmt: skip

    def fake_happy(job: HappyJob, toolbox: Toolbox, manifest: Any, step: str) -> None:
        job.prefix.parent.mkdir(parents=True, exist_ok=True)
        _write_vcf(Path(f"{job.prefix}.vcf.gz"), ["TRUTH", "QUERY"], happy_rows)
        Path(f"{job.prefix}.extended.csv").write_text(EXTENDED)
        assert "planted" in job.strata.read_text()

    return run, pipeline, {"happy": fake_happy}


def test_compare_sample_end_to_end(tmp_path: Path, reference: Reference) -> None:
    run, pipeline, kw = _run_dir(tmp_path, reference)
    args = {"reference": reference, "strata": [], "toolbox": Toolbox(), "manifest": None,
            "params": CompareParams()}  # fmt: skip
    result = compare_sample(run, pipeline, "s", **args, **kw)
    classes = {r.variant_id: r.outcome_class for r in result.records}
    assert classes == {
        "v1_plp": OutcomeClass.CORRECT,
        "v2_plp_del": OutcomeClass.NOT_CALLED,
        "v3_blb": OutcomeClass.MISCLASSIFIED,
        "v4_calling": OutcomeClass.CORRECT,
        "v5_delins": OutcomeClass.NOT_ANNOTATED,
    }
    fr = {(r.chrom, r.pos): r.evidence for r in result.false_reports}
    assert fr[("chr1", 2001)]["planted_negative"] is True
    assert fr[("chr1", 2001)]["planted_variant_id"] == "v3_blb"
    assert fr[("chr1", 4001)]["call"] == "true_variant"
    assert fr[("chr1", 7001)]["call"] == "fp_call"
    assert fr[("chrX", 501)]["call"] == "outside_scored_regions"
    assert fr[("chr1", 9301)]["call"] == "no_exact_match"
    assert ("chr1", 9001) not in fr  # linked to v5, not counted again
    split = {r.variant_id: r.evidence.get("split_report") for r in result.records}
    assert split["v5_delins"].startswith("chr1:9001:") and split["v5_delins"].endswith("T")
    assert split["v1_plp"] is None
    held = {r.variant_id: r.evidence["held_out"] for r in result.records}
    assert held == {
        "v1_plp": True, "v2_plp_del": None, "v3_blb": False, "v4_calling": None, "v5_delins": None
    }  # fmt: skip
    assert len(fr) == 5
    out = compare_dir(run, pipeline)
    calling = json.loads((out / "s.calling.json").read_text())
    assert {(m["type"], m["stratum"]) for m in calling} == {
        ("SNP", "*"),
        ("INDEL", "*"),
        ("SNP", "planted"),
    }
    first = (out / "s.outcomes.json").read_bytes(), write_outcomes_tsv(out).read_bytes()
    compare_sample(run, pipeline, "s", **args, **kw)
    assert ((out / "s.outcomes.json").read_bytes(), write_outcomes_tsv(out).read_bytes()) == first
    rows = write_outcomes_tsv(out).read_text().splitlines()
    assert len(rows) == 1 + 5 + 5 and rows[0].startswith("pipeline\tsample_id\tvariant_id")


def test_compare_sample_records_manifest(tmp_path: Path, reference: Reference) -> None:
    run, pipeline, kw = _run_dir(tmp_path, reference)
    manifest = ManifestBuilder(run, config_path=tmp_path / "p.yaml", config_sha256="0" * 64, seed=1)
    compare_sample(run, pipeline, "s", reference=reference, strata=[], toolbox=Toolbox(),
                   manifest=manifest, params=CompareParams(), **kw)  # fmt: skip
    roles = {o.role for o in manifest.manifest.outputs}
    assert {"outcomes:p1:s", "calling:p1:s"} <= roles


def test_happy_gender_argv(tmp_path: Path) -> None:
    """compare passes --gender to hap.py only when asked."""
    job = HappyJob(
        truth=tmp_path / "t.vcf.gz", query=tmp_path / "q.vcf.gz", reference=tmp_path / "r.fa",
        regions=tmp_path / "r.bed", strata=tmp_path / "s.tsv", prefix=tmp_path / "p", threads=2,
    )  # fmt: skip
    assert "--gender" not in happy_args(job)
    male = HappyJob(**{**job.__dict__, "gender": Sex.MALE})
    assert happy_args(male)[-2:] == ["--gender", "male"]
