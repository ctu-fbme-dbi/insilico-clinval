"""ClinGen eRepo: parsing, held_out, selection and pairing."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from insilico_clinval.config import SuiteError, load_suite
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.enums import Classification
from insilico_clinval.schemas.selection import RuleSummary, SelectionResult
from insilico_clinval.truth import TruthError
from insilico_clinval.truth.genotype import assign_genotypes
from insilico_clinval.variants.erepo import (
    ErepoError,
    ErepoMismatchError,
    HeldOutReference,
    load_erepo,
    load_held_out_reference,
    reference_ids,
)
from insilico_clinval.variants.report import selected_tsv, selection_report
from insilico_clinval.variants.select import (
    SelectionError,
    select_variants,
    selection_failures,
    variants_per_sample,
)
from tests.conftest import EREPO_SUITE, FIXTURES, rule

REFERENCE_UUIDS = {
    "EREPO-0001",
    "EREPO-0004",
    "EREPO-0005",
    "EREPO-0014",
    "EREPO-0016",
    "EREPO-0020",
    "EREPO-0023",
    "EREPO-0030",
    "EREPO-0032",
}


@pytest.fixture
def result(reference: Reference) -> SelectionResult:
    return select_variants(load_suite(EREPO_SUITE), reference)


def _select(
    data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path], reference: Reference
) -> SelectionResult:
    return select_variants(load_suite(write_suite(data)), reference)


def _summary(result: SelectionResult, rule_id: str) -> RuleSummary:
    return next(r for r in result.rules if r.rule_id == rule_id)


def test_load_erepo_merges_skips_and_normalises(reference: Reference) -> None:
    ids = reference_ids(FIXTURES / "erepo_v2.2.0_mini.tsv")
    assert ids >= REFERENCE_UUIDS and "CA900001" in ids
    loaded = load_erepo(
        FIXTURES / "erepo_mini.tsv",
        release="r1",
        reference=reference,
        held_out_reference=HeldOutReference("v2.2.0", ids=frozenset(ids)),
    )
    assert loaded.rows == 27
    assert loaded.skipped == {
        "discordant_classifications": 1,
        "no_grch38_hgvs": 1,
        "retracted": 1,
        "unsupported_hgvs": 1,
    }
    by_uuid = {a.primary.uuid: a for a in loaded.alleles}
    merged = by_uuid["EREPO-0005"]
    assert merged.uuids == ("EREPO-0005", "EREPO-0006")
    assert merged.classification == Classification.PATHOGENIC
    assert [c.id for c in merged.conditions] == ["MONDO:0000101", "MONDO:0000111"]
    assert not merged.held_out and merged.release_first_seen == "v2.2.0"
    # HGVS g.9505del is 3'-shifted in the chr2:9500-9505 homopolymer; VCF is left-aligned
    homopolymer = by_uuid["EREPO-0021"]
    assert (homopolymer.chrom, homopolymer.pos, homopolymer.ref, homopolymer.alt) == (
        "chr2",
        9499,
        "CA",
        "C",
    )
    assert homopolymer.held_out and homopolymer.release_first_seen == "r1"
    assert "EREPO-0007" not in by_uuid and "EREPO-0008" not in by_uuid  # discordant P vs B
    assert by_uuid["EREPO-0001"].primary.met_criteria == ("PVS1", "PM2_Supporting")


def test_held_out_reference_vcf(tmp_path: Path, reference: Reference) -> None:
    """A VCF reference matches by left-normalized allele (here a right-shifted deletion)."""
    seq = reference.fetch("chr2", 9497, 9506)  # chr2:9498-9506, homopolymer at 9500-9505
    vcf = tmp_path / "benchmark.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr2,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr2\t9504\t.\t{seq[6:8]}\t{seq[6]}\t.\tPASS\t.\n"
    )
    ref = load_held_out_reference(vcf, "v2.2.0", reference)
    assert ref.alleles == {("chr2", 9499, "CA", "C")} and not ref.ids
    # an unanchored deletion with an empty ALT (BIAS-2015 benchmark style) is the same allele
    unanchored = tmp_path / "unanchored.vcf"
    unanchored.write_text(vcf.read_text().replace(f"\t{seq[6:8]}\t{seq[6]}\t", f"\t{seq[7]}\t\t"))
    unanchored.write_text(unanchored.read_text().replace("\t9504\t", "\t9505\t"))
    assert load_held_out_reference(unanchored, "v2.2.0", reference).alleles == ref.alleles
    # a mangled indel (wrong length) near the allele still counts as present; an SNV does not
    mangled = HeldOutReference("v2.2.0", indel_sites=frozenset({("chr2", 9506)}))
    assert mangled.contains(("chr2", 9499, "CAAAA", "C"), [])
    assert not mangled.contains(("chr2", 9499, "CA", "C"), [])  # span 9499-9500 + 3 bp < 9506
    assert not mangled.contains(("chr2", 9506, "A", "G"), [])
    loaded = load_erepo(
        FIXTURES / "erepo_mini.tsv", release="r1", reference=reference, held_out_reference=ref
    )
    by_uuid = {a.primary.uuid: a for a in loaded.alleles}
    assert not by_uuid["EREPO-0021"].held_out
    assert by_uuid["EREPO-0021"].release_first_seen == "v2.2.0"
    assert by_uuid["EREPO-0001"].held_out
    vcf.write_text(vcf.read_text().replace(f"\t{seq[6:8]}\t", "\tGG\t"))
    with pytest.raises(ErepoError, match="do not match the reference FASTA"):
        load_held_out_reference(vcf, "v2.2.0", reference)
    tsv = load_held_out_reference(FIXTURES / "erepo_v2.2.0_mini.tsv", "v2.2.0", reference)
    assert tsv.ids >= REFERENCE_UUIDS and not tsv.alleles


def test_missing_column_and_unknown_assertion(tmp_path: Path, reference: Reference) -> None:
    src = (FIXTURES / "erepo_mini.tsv").read_text().splitlines()
    bad_header = tmp_path / "no_uuid.tsv"
    bad_header.write_text("\n".join([src[0].replace("\tUuid", "\tUUID"), *src[1:]]) + "\n")
    with pytest.raises(ErepoError, match=r"lacks columns \['Uuid'\]"):
        load_erepo(bad_header, release="r", reference=reference, held_out_reference=None)
    bad_value = tmp_path / "assertion.tsv"
    bad_value.write_text("\n".join([src[0], src[1].replace("\tPathogenic\t", "\tHarmful\t")]))
    with pytest.raises(ErepoError, match="unknown Assertion 'Harmful'"):
        load_erepo(bad_value, release="r", reference=reference, held_out_reference=None)
    no_ids = tmp_path / "ref.tsv"
    no_ids.write_text("Gene\tAssertion\nGENEA\tBenign\n")
    with pytest.raises(ErepoError, match="needs a 'Uuid' or 'Allele Registry Id' column"):
        reference_ids(no_ids)


def test_stated_reference_mismatch_is_an_error(reference: Reference) -> None:
    """ERepo coordinates: a stated base that differs from the FASTA."""
    with pytest.raises(ErepoMismatchError, match=r"EREPO-0099.*stated reference A .*\(C\)"):
        load_erepo(
            FIXTURES / "erepo_badref.tsv", release="r", reference=reference, held_out_reference=None
        )


def test_erepo_provenance_and_held_out(result: SelectionResult) -> None:
    """Source, release, record id, VCEP and held_out on every eRepo variant."""
    erepo = [v for v in result.variants if v.source.database == "ClinGen eRepo"]
    assert len(erepo) == 12
    for v in erepo:
        assert v.clingen is not None
        assert v.source.release == "fixture-2026-09"
        assert v.source.record_id == v.clingen.erepo_id
        assert v.clingen.vcep == f"Fixture Gene {v.gene[-1]} VCEP"  # type: ignore[index]
        in_reference = v.clingen.erepo_id in REFERENCE_UUIDS
        assert v.held_out is not in_reference
        assert v.clingen.erepo_release_first_seen == (
            "v2.2.0" if in_reference else "fixture-2026-09"
        )
    assert {v.held_out for v in erepo} == {True, False}
    synthetic = next(v for v in result.variants if v.category.value == "synthetic_novel")
    assert synthetic.clingen is None and synthetic.held_out is False


def test_same_genes_and_allow_list(result: SelectionResult) -> None:
    """ERepo: VUS and B/LB never pass a P/LP list; P/LP and B/LB use the same genes."""
    plp = _summary(result, "erepo_plp")
    blb = _summary(result, "erepo_blb")
    assert plp.satisfied and blb.satisfied
    assert plp.selected == blb.selected == {"GENEA": 2, "GENEB": 2, "GENEX": 2}
    assert all(v.source.record_id != "EREPO-0013" for v in result.variants)  # VUS
    counts = dict(plp.filter_counts)
    assert counts["classification_allow"] < counts["max_indel_length"] < counts["variant_type"]
    classes = {v.clingen.classification for v in result.variants if v.rule_id == "erepo_plp"}  # type: ignore[union-attr]
    assert classes <= {Classification.PATHOGENIC, Classification.LIKELY_PATHOGENIC}
    assert all(v.source.record_id != "EREPO-0012" for v in result.variants)  # 25 bp > 20


def test_pairing_report_per_gene(result: SelectionResult) -> None:
    """Per-gene positive / negative counts with held-out split and VCEP."""
    report = selection_report(load_suite(EREPO_SUITE), result)
    assert "- ClinGen eRepo release: fixture-2026-09 (held-out reference v2.2.0)" in report
    assert "- background: GIAB_HG002 fixture-v4.2.1 (male) in samples positive, negative" in report
    rows = {
        line.split("|")[1].strip(): [c.strip() for c in line.split("|")[2:-1]]
        for line in report.splitlines()
        if line.startswith("| GENE")
    }
    assert set(rows) == {"GENEA", "GENEB", "GENEX"}
    per_sample = variants_per_sample(result)
    for gene, (pos, pos_ho, neg, neg_ho, vcep, note) in rows.items():
        assert (pos, neg, note) == ("2", "2", "")
        assert vcep == f"Fixture Gene {gene[-1]} VCEP"
        assert int(pos_ho) == sum(v.held_out for v in per_sample["positive"] if v.gene == gene)
        assert int(neg_ho) == sum(v.held_out for v in per_sample["negative"] if v.gene == gene)
    assert "below the 59 per category recommended by Roy et al. 2018" in report
    tsv = selected_tsv(result).splitlines()
    header = tsv[0].split("\t")
    assert {"vcep", "erepo_classification", "met_criteria", "held_out"} <= set(header)


def test_asymmetry_is_flagged(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    rule(erepo_suite_data, "erepo_blb")["genes"] = ["GENEA", "GENEB"]
    suite = load_suite(write_suite(erepo_suite_data))
    res = select_variants(suite, reference)
    assert "| GENEX | 2 | " in selection_report(suite, res)
    assert "**asymmetric**" in selection_report(suite, res)


def test_too_close_erepo_candidate_rejected(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    """EREPO-0004 (chr1:400) is 99 bp from EREPO-0001 (chr1:300) -> rejected, logged."""
    r = rule(erepo_suite_data, "erepo_plp")
    for key in ("genes", "per_gene"):
        r.pop(key)
    r["erepo_ids"] = ["EREPO-0001", "EREPO-0004", "EREPO-0006"]
    summary = _summary(_select(erepo_suite_data, write_suite, reference), "erepo_plp")
    assert summary.selected == {"eRepo EREPO-0001": 1, "eRepo EREPO-0004": 0, "eRepo EREPO-0006": 1}
    (rejection,) = summary.rejections
    assert rejection.variant_id.startswith("curated:chr1:400:")
    assert "99 bp from curated:chr1:300:" in rejection.reason
    assert "min_distance_bp=300" in rejection.reason
    assert any(m.startswith("eRepo EREPO-0004: requested 1, selected 0") for m in summary.missing)


def test_held_out_filters(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    r = rule(erepo_suite_data, "erepo_blb")
    r.update(held_out="only", required=False)
    res = _select(erepo_suite_data, write_suite, reference)
    blb = [v for v in res.variants if v.rule_id == "erepo_blb"]
    assert blb and all(v.held_out for v in blb)
    r["held_out"] = "exclude"
    res = _select(erepo_suite_data, write_suite, reference)
    blb = [v for v in res.variants if v.rule_id == "erepo_blb"]
    assert blb and not any(v.held_out for v in blb)
    assert not _summary(res, "erepo_blb").satisfied  # GENEA/GENEB/GENEX have one B/LB each


def test_held_out_rule_needs_reference(
    erepo_suite_data: dict[str, Any], write_suite: Callable[[dict[str, Any]], Path]
) -> None:
    erepo_suite_data["sources"]["clingen_erepo"].pop("held_out_reference")
    rule(erepo_suite_data, "erepo_plp")["held_out"] = "only"
    with pytest.raises(SuiteError, match=r"held_out: only needs sources\.clingen_erepo"):
        load_suite(write_suite(erepo_suite_data))


def test_without_reference_held_out_is_false_with_warning(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    erepo_suite_data["sources"]["clingen_erepo"].pop("held_out_reference")
    res = _select(erepo_suite_data, write_suite, reference)
    assert not any(v.held_out for v in res.variants)
    assert any("no held_out_reference configured" in w for w in res.warnings)


def test_unsatisfied_required_erepo_rule(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    """A required eRepo rule with a missing record id makes the selection fail."""
    r = rule(erepo_suite_data, "erepo_blb")
    for key in ("genes", "per_gene"):
        r.pop(key)
    r["erepo_ids"] = ["EREPO-0014", "EREPO-9999"]
    res = _select(erepo_suite_data, write_suite, reference)
    failures = selection_failures(res, ["positive", "negative"])
    assert failures == [
        "required rule erepo_blb is unsatisfied: eRepo EREPO-9999: requested 1, selected 0 "
        "(0 candidates passed all filters, 0 rejected by deduplication)"
    ]


def test_erepo_indel_ladder_shortfall(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    """ERepo: the fixture has only 1-5 bp indels; longer bins report the shortfall."""
    erepo_suite_data["selection"].append(
        {
            "id": "erepo_ladder",
            "category": "difficult",
            "source": "clingen_erepo",
            "classification_allow": ["pathogenic", "likely_pathogenic", "benign", "likely_benign"],
            "ladder_per_bin": 1,
            "expected": {"reportable": False, "level": "calling"},
            "samples": ["positive"],
        }
    )
    summary = _summary(_select(erepo_suite_data, write_suite, reference), "erepo_ladder")
    assert summary.requested == {"1-5bp": 1, "6-10bp": 1, "11-20bp": 1}
    assert summary.selected["6-10bp"] == summary.selected["11-20bp"] == 0
    assert any(m.startswith("6-10bp: requested 1, selected 0") for m in summary.missing)


def test_compound_het_on_erepo_rule(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
) -> None:
    """Compound heterozygotes (trans) also work for eRepo rules (female, no background)."""
    erepo_suite_data.pop("background")
    for s in erepo_suite_data["samples"]:
        s["sex"] = "female"
    r = rule(erepo_suite_data, "erepo_plp")
    r["compound_het"] = {"phase": "trans"}
    suite = load_suite(write_suite(erepo_suite_data))
    res = select_variants(suite, reference)
    sample = next(s for s in suite.config.samples if s.sample_id == "positive")
    genotypes, _ = assign_genotypes(suite.config, sample, variants_per_sample(res)["positive"])
    pairs: dict[int, set[tuple[int, ...]]] = {}
    for g in genotypes:
        if g.variant.rule_id == "erepo_plp":
            assert g.phase_set is not None
            pairs.setdefault(g.phase_set, set()).add(g.gt)
    assert len(pairs) == 3 and all(gts == {(1, 0), (0, 1)} for gts in pairs.values())


def test_background_sex_and_giab_conflict(
    erepo_suite_data: dict[str, Any],
    suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
) -> None:
    erepo_suite_data["samples"][0]["sex"] = "female"
    with pytest.raises(SuiteError, match="sex female differs from the GIAB_HG002 background"):
        load_suite(write_suite(erepo_suite_data))
    # a GIAB rule drawing from the very VCF used as background of its samples
    erepo_suite_data["samples"][0]["sex"] = "male"
    erepo_suite_data["background"]["vcf"] = suite_data["sources"]["giab"]["truth_vcf"]
    erepo_suite_data["sources"]["giab"] = suite_data["sources"]["giab"]
    erepo_suite_data["selection"].append(rule(suite_data, "standard_giab"))
    with pytest.raises(SuiteError, match="rule standard_giab: selects from the background VCF"):
        load_suite(write_suite(erepo_suite_data))
    erepo_suite_data["samples"][0]["background"] = False
    erepo_suite_data["samples"][1]["background"] = False
    load_suite(write_suite(erepo_suite_data))  # allowed without background


def test_selection_error_types(
    erepo_suite_data: dict[str, Any],
    write_suite: Callable[[dict[str, Any]], Path],
    reference: Reference,
    tmp_path: Path,
) -> None:
    broken = tmp_path / "broken.tsv"
    broken.write_text("Uuid\tAssertion\nX\tBenign\n")
    erepo_suite_data["sources"]["clingen_erepo"]["file"] = str(broken)
    with pytest.raises(SelectionError, match="eRepo export lacks columns"):
        _select(erepo_suite_data, write_suite, reference)
    assert issubclass(ErepoMismatchError, ValueError)
    assert not issubclass(ErepoMismatchError, TruthError)  # mapped to exit 3 by the CLI
