"""Zygosity, haplotype and phase-set assignment for one sample.

Rules (docs/schemas.md, SampleSpec):

- zygosity: sample override > rule ``zygosity`` > sample ``default_zygosity``;
- male samples: chrX/chrY outside the pseudo-autosomal regions are haploid -> ``hemi``, GT ``1``;
  a chrY variant in a female sample is an error;
- het haplotype: ``derive_choice(seed, "haplotype", f"{sample}|{variant_id}", 2)``
  (0 -> haplotype 1, GT ``1|0``; 1 -> haplotype 2, GT ``0|1``);
- phase sets (compound heterozygotes): from rule ``compound_het`` (consecutive pairs by position
  within a gene) and sample ``phase_sets``; ``PS`` = position of the first member; the first
  member's haplotype is ``derive_choice(seed, "phase_set", f"{sample}|{chrom}:{PS}", 2)``, cis
  members share it, the trans partner takes the other one.
"""

from __future__ import annotations

from dataclasses import dataclass

from insilico_clinval.constants import GRCH38_PAR_REGIONS
from insilico_clinval.rng import derive_choice
from insilico_clinval.schemas.enums import Haplotype, Phase, Sex, Zygosity
from insilico_clinval.schemas.sample import SampleVariant
from insilico_clinval.schemas.suite import (
    ClinGenRule,
    ClinVarRule,
    Region,
    SampleConfig,
    SuiteConfig,
)
from insilico_clinval.schemas.variant import TestVariant
from insilico_clinval.truth import TruthError


@dataclass(frozen=True)
class Genotype:
    variant: TestVariant
    zygosity: Zygosity
    haplotype: Haplotype
    phase_set: int | None
    notes: tuple[str, ...] = ()

    @property
    def gt(self) -> tuple[int, ...]:
        if self.zygosity == Zygosity.HEMI:
            return (1,)
        if self.zygosity == Zygosity.HOM:
            return (1, 1)
        return (1, 0) if self.haplotype == Haplotype.ONE else (0, 1)

    def to_sample_variant(self) -> SampleVariant:
        return SampleVariant(
            variant_id=self.variant.variant_id,
            zygosity=self.zygosity,
            haplotype=self.haplotype,
            phase_set=self.phase_set,
        )


def par_regions(cfg: SuiteConfig) -> tuple[Region, ...]:
    if cfg.reference.par_regions is not None:
        return tuple(cfg.reference.par_regions)
    return GRCH38_PAR_REGIONS


def is_haploid(chrom: str, pos: int, sex: Sex, pars: tuple[Region, ...]) -> bool:
    if sex == Sex.FEMALE or chrom not in ("chrX", "chrY"):
        return False
    return not any(r.chrom == chrom and r.start <= pos <= r.end for r in pars)


def _hap(index: int) -> Haplotype:
    return Haplotype.ONE if index == 0 else Haplotype.TWO


def _other(h: Haplotype) -> Haplotype:
    return Haplotype.TWO if h == Haplotype.ONE else Haplotype.ONE


def assign_genotypes(
    cfg: SuiteConfig, sample: SampleConfig, variants: list[TestVariant]
) -> tuple[list[Genotype], list[str]]:
    """Return the genotypes of ``sample`` (input order) and warnings."""
    seed = cfg.suite.seed
    pars = par_regions(cfg)
    rules = {r.id: r for r in cfg.selection}
    by_id = {v.variant_id: v for v in variants}
    overrides = {o.variant_id: o.zygosity for o in sample.overrides}
    warnings: list[str] = []
    errors: list[str] = []

    unknown = sorted(set(overrides) - set(by_id))
    if unknown:
        errors.append(f"overrides reference variants not in sample {sample.sample_id}: {unknown}")

    zygosity: dict[str, Zygosity] = {}
    notes: dict[str, list[str]] = {v.variant_id: [] for v in variants}
    for v in variants:
        z = overrides.get(v.variant_id) or rules[v.rule_id].zygosity or sample.default_zygosity
        if sample.sex == Sex.FEMALE and v.chrom == "chrY":
            errors.append(f"{v.variant_id}: chrY variant in female sample {sample.sample_id}")
            continue
        if is_haploid(v.chrom, v.pos, sample.sex, pars):
            if z != Zygosity.HEMI:
                notes[v.variant_id].append(f"{z.value} -> hemi: male sample, {v.chrom} outside PAR")
            z = Zygosity.HEMI
        elif z == Zygosity.HEMI:
            errors.append(f"{v.variant_id}: hemi requested on a diploid locus")
            continue
        zygosity[v.variant_id] = z

    # Phase sets: lists of variant ids with their phase.
    groups: list[tuple[list[str], Phase, str]] = []
    for rule in cfg.selection:
        if not isinstance(rule, (ClinGenRule, ClinVarRule)) or rule.compound_het is None:
            continue
        if sample.sample_id not in rule.samples:
            continue
        by_gene: dict[str, list[TestVariant]] = {}
        for v in variants:
            if v.rule_id == rule.id and v.gene is not None:
                by_gene.setdefault(v.gene, []).append(v)
        for gene in sorted(by_gene):
            members = sorted(by_gene[gene], key=lambda x: x.pos)
            for i in range(0, len(members) - 1, 2):
                groups.append(
                    (
                        [members[i].variant_id, members[i + 1].variant_id],
                        rule.compound_het.phase,
                        f"rule {rule.id}, gene {gene}",
                    )
                )
            if len(members) % 2:
                warnings.append(
                    f"sample {sample.sample_id}: rule {rule.id} gene {gene} has an odd number of "
                    f"variants; {members[-1].variant_id} is not in a compound-het pair"
                )
    for spec in sample.phase_sets:
        missing = [vid for vid in spec.variants if vid not in by_id]
        if missing:
            errors.append(
                f"phase_sets reference variants not in sample {sample.sample_id}: {missing}"
            )
            continue
        groups.append((list(spec.variants), spec.phase, "sample phase_sets"))

    phase_of: dict[str, tuple[int, Haplotype, str]] = {}
    for ids, phase, origin in groups:
        members = sorted((by_id[i] for i in ids), key=lambda x: x.pos)
        chroms = {m.chrom for m in members}
        problems = [
            f"{m.variant_id} is {zygosity[m.variant_id].value}"
            for m in members
            if m.variant_id in zygosity and zygosity[m.variant_id] != Zygosity.HET
        ]
        if len(chroms) != 1:
            problems.append("members on different chromosomes")
        if any(m.variant_id in phase_of for m in members):
            problems.append("a member is already in another phase set")
        if problems:
            errors.append(f"phase set ({origin}) {ids}: " + "; ".join(problems))
            continue
        ps = members[0].pos
        first = _hap(
            derive_choice(seed, "phase_set", f"{sample.sample_id}|{members[0].chrom}:{ps}", 2)
        )
        for k, m in enumerate(members):
            h = first if (phase == Phase.CIS or k == 0) else _other(first)
            phase_of[m.variant_id] = (ps, h, f"compound het {phase.value} ({origin})")

    if errors:
        raise TruthError(f"sample {sample.sample_id}:\n  " + "\n  ".join(errors))

    genotypes: list[Genotype] = []
    for v in variants:
        z = zygosity[v.variant_id]
        phase_set: int | None = None
        if z == Zygosity.HEMI:
            h = Haplotype.ONE
        elif z == Zygosity.HOM:
            h = Haplotype.BOTH
        elif v.variant_id in phase_of:
            phase_set, h, note = phase_of[v.variant_id]
            notes[v.variant_id].append(note)
        else:
            h = _hap(derive_choice(seed, "haplotype", f"{sample.sample_id}|{v.variant_id}", 2))
        genotypes.append(
            Genotype(
                variant=v,
                zygosity=z,
                haplotype=h,
                phase_set=phase_set,
                notes=tuple(notes[v.variant_id]),
            )
        )
    return genotypes, warnings
