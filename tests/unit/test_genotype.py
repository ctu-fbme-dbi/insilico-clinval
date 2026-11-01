"""Zygosity, haplotype and phase-set assignment."""

from __future__ import annotations

from typing import Any

import pytest

from insilico_clinval.schemas.enums import Haplotype, VariantCategory, VariantType, Zygosity
from insilico_clinval.schemas.suite import SuiteConfig
from insilico_clinval.schemas.variant import Source, TestVariant
from insilico_clinval.truth import TruthError
from insilico_clinval.truth.genotype import assign_genotypes, is_haploid, par_regions
from tests.conftest import rule


def _tv(chrom: str, pos: int, rule_id: str, gene: str | None = "GENEA") -> TestVariant:
    return TestVariant(
        variant_id=f"curated:{chrom}:{pos}:A:G",
        rule_id=rule_id,
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="G",
        variant_type=VariantType.SNV,
        indel_length=0,
        category=VariantCategory.CURATED,
        gene=gene,
        source=Source(database="ClinVar", release="x", record_id=str(pos)),
        clinvar=None,
        strata=[],
        selection_reason="test",
    )


def _cfg(data: dict[str, Any]) -> SuiteConfig:
    return SuiteConfig.model_validate(data)


def _gt(genotypes: list[Any], pos: int) -> Any:
    return next(g for g in genotypes if g.variant.pos == pos)


def test_het_is_phased_and_deterministic(suite_data: dict[str, Any]) -> None:
    rule(suite_data, "clinvar_plp_panel").pop("compound_het")
    cfg = _cfg(suite_data)
    variants = [_tv("chr1", p, "clinvar_plp_panel", gene=f"G{p}") for p in (100, 2000, 4000, 6000)]
    g1, _ = assign_genotypes(cfg, cfg.samples[0], variants)
    g2, _ = assign_genotypes(cfg, cfg.samples[0], variants)
    assert [g.gt for g in g1] == [g.gt for g in g2]
    assert all(g.gt in ((1, 0), (0, 1)) and g.zygosity == Zygosity.HET for g in g1)
    assert all(g.phase_set is None for g in g1)


def test_hom_from_rule_and_override(suite_data: dict[str, Any]) -> None:
    rule(suite_data, "clinvar_plp_panel").pop("compound_het")
    rule(suite_data, "clinvar_plp_panel")["zygosity"] = "hom"
    suite_data["samples"][0]["overrides"] = [
        {"variant_id": "curated:chr1:3000:A:G", "zygosity": "het"}
    ]
    cfg = _cfg(suite_data)
    variants = [
        _tv("chr1", 1000, "clinvar_plp_panel", "G1"),
        _tv("chr1", 3000, "clinvar_plp_panel", "G2"),
    ]
    genotypes, _ = assign_genotypes(cfg, cfg.samples[0], variants)
    assert _gt(genotypes, 1000).gt == (1, 1)
    assert _gt(genotypes, 1000).haplotype == Haplotype.BOTH
    assert _gt(genotypes, 3000).zygosity == Zygosity.HET


def test_compound_het_trans_from_rule(suite_data: dict[str, Any]) -> None:
    cfg = _cfg(suite_data)
    variants = [_tv("chr1", 1000, "clinvar_plp_panel"), _tv("chr1", 3000, "clinvar_plp_panel")]
    genotypes, warnings = assign_genotypes(cfg, cfg.samples[0], variants)
    a, b = _gt(genotypes, 1000), _gt(genotypes, 3000)
    assert a.phase_set == b.phase_set == 1000
    assert {a.gt, b.gt} == {(1, 0), (0, 1)}
    assert not warnings


def test_compound_het_cis_from_sample(suite_data: dict[str, Any]) -> None:
    rule(suite_data, "clinvar_plp_panel").pop("compound_het")
    suite_data["samples"][0]["phase_sets"] = [
        {"variants": ["curated:chr1:3000:A:G", "curated:chr1:1000:A:G"], "phase": "cis"}
    ]
    cfg = _cfg(suite_data)
    variants = [_tv("chr1", 1000, "clinvar_plp_panel"), _tv("chr1", 3000, "clinvar_plp_panel")]
    genotypes, _ = assign_genotypes(cfg, cfg.samples[0], variants)
    a, b = _gt(genotypes, 1000), _gt(genotypes, 3000)
    assert a.phase_set == b.phase_set == 1000
    assert a.gt == b.gt


def test_odd_compound_het_member_warns(suite_data: dict[str, Any]) -> None:
    cfg = _cfg(suite_data)
    variants = [_tv("chr1", p, "clinvar_plp_panel") for p in (1000, 3000, 5000)]
    genotypes, warnings = assign_genotypes(cfg, cfg.samples[0], variants)
    assert _gt(genotypes, 5000).phase_set is None
    assert any("odd number" in w for w in warnings)


def test_male_chrx_is_hemizygous_outside_par(suite_data: dict[str, Any]) -> None:
    cfg = _cfg(suite_data)
    male = cfg.samples[1]
    variants = [
        _tv("chrX", 4000, "clinvar_blb_panel", "G1"),
        _tv("chrX", 500, "clinvar_blb_panel", "G2"),
    ]
    genotypes, _ = assign_genotypes(cfg, male, variants)
    hemi, par = _gt(genotypes, 4000), _gt(genotypes, 500)
    assert hemi.zygosity == Zygosity.HEMI and hemi.gt == (1,)
    assert "hemi" in hemi.notes[0]
    assert par.zygosity == Zygosity.HET and len(par.gt) == 2


def test_chry_in_female_is_an_error(suite_data: dict[str, Any]) -> None:
    cfg = _cfg(suite_data)
    with pytest.raises(TruthError, match="chrY"):
        assign_genotypes(cfg, cfg.samples[0], [_tv("chrY", 1500, "clinvar_plp_panel", "G1")])


def test_compound_het_on_hemizygous_locus_is_an_error(suite_data: dict[str, Any]) -> None:
    rule(suite_data, "clinvar_blb_panel")["compound_het"] = {"phase": "trans"}
    cfg = _cfg(suite_data)
    variants = [_tv("chrX", 2000, "clinvar_blb_panel"), _tv("chrX", 4000, "clinvar_blb_panel")]
    with pytest.raises(TruthError, match="hemi"):
        assign_genotypes(cfg, cfg.samples[1], variants)


def test_default_pars_are_grch38(suite_data: dict[str, Any]) -> None:
    suite_data["reference"].pop("par_regions")
    pars = par_regions(_cfg(suite_data))
    from insilico_clinval.schemas.enums import Sex

    assert not is_haploid("chrX", 2_000_000, Sex.MALE, pars)  # PAR1
    assert is_haploid("chrX", 50_000_000, Sex.MALE, pars)
    assert not is_haploid("chrX", 50_000_000, Sex.FEMALE, pars)
