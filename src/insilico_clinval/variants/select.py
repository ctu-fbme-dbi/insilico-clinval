"""Candidate selection: filters, ranking, per-bucket pools and suite-wide deduplication.

Algorithm (documented in docs/architecture.md):

1. Each source is streamed once; every record is evaluated against every rule of that source.
   Per rule, a chain of filters is applied and the number of candidates remaining after each
   filter is counted (reported in ``selection_report.md``).
2. Passing candidates go to a *bucket* of the rule (gene for ``per_gene``, ladder bin, variant
   type for GIAB ``counts``, requested id, or ``total``). Each bucket keeps a bounded pool of the
   best-ranked candidates.
3. Ranking: ClinVar review stars descending, then ``derive_int(seed, "rank", variant_id)``
   ascending (GIAB: hash only; explicit variants: listed order).
4. Rules are processed in YAML order, buckets in configuration order; candidates are accepted
   greedily in rank order. A candidate whose REF span overlaps an accepted variant or lies
   closer than ``min_distance_bp`` (gap between REF spans) is rejected with a logged reason and
   the next ranked candidate is tried.
"""

from __future__ import annotations

import heapq
import logging
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field, replace
from pathlib import Path

import pysam

from insilico_clinval._exec import elapsed
from insilico_clinval.config import LoadedSuite
from insilico_clinval.genome.alleles import (
    INDEL_TYPES,
    classify,
    event_length,
    indel_length,
    is_simple_allele,
    ref_span0,
)
from insilico_clinval.genome.bed import IntervalSet
from insilico_clinval.genome.reference import Reference
from insilico_clinval.rng import derive_int
from insilico_clinval.schemas import SCHEMA_VERSION
from insilico_clinval.schemas.enums import VariantCategory, VariantType
from insilico_clinval.schemas.selection import Rejection, RuleSummary, SelectionResult
from insilico_clinval.schemas.suite import ClinGenRule, ClinVarRule, ExplicitRule, GiabRule, Rule
from insilico_clinval.schemas.variant import ClinGenInfo, ClinVarInfo, Source, TestVariant
from insilico_clinval.variants.clinvar import ClinVarAnnotation, parse_record, passes_allow_list
from insilico_clinval.variants.erepo import DATABASE as EREPO_DATABASE
from insilico_clinval.variants.erepo import (
    ErepoAllele,
    ErepoError,
    HeldOutReference,
    load_erepo,
    load_held_out_reference,
)

log = logging.getLogger(__name__)
# Phase messages only (shown on the console); per-candidate details go to ``log``.
progress = logging.getLogger(f"{__name__}.progress")

# Minimum slack of candidates kept per bucket beyond the requested count.
POOL_SLACK = 1000


class SelectionError(RuntimeError):
    """Selection could not be carried out (invalid input or an exhausted candidate pool)."""


@dataclass(frozen=True)
class Candidate:
    chrom: str
    pos: int
    ref: str
    alt: str
    variant_type: VariantType
    event_length: int
    gene: str | None
    source: Source
    clinvar: ClinVarAnnotation | None
    erepo: ErepoAllele | None = None
    note: str = ""

    @property
    def stars(self) -> int:
        return self.clinvar.stars if self.clinvar is not None else 0

    def variant_id(self, category: VariantCategory) -> str:
        return TestVariant.make_id(category, self.chrom, self.pos, self.ref, self.alt)


class Pool:
    """The ``capacity`` best-ranked candidates of one bucket (min key = best)."""

    def __init__(self, capacity: int):
        self.capacity = capacity
        self.truncated = False
        self.seen = 0
        self._heap: list[tuple[tuple[int, ...], int, Candidate]] = []

    def push(self, key: tuple[int, ...], candidate: Candidate) -> None:
        self.seen += 1
        # Max-heap on key through negation; ``seen`` breaks exact ties by input order.
        item = (tuple(-k for k in key), -self.seen, candidate)
        if len(self._heap) < self.capacity:
            heapq.heappush(self._heap, item)
        else:
            self.truncated = True
            heapq.heappushpop(self._heap, item)

    def ranked(self) -> list[Candidate]:
        return [c for _, _, c in sorted(self._heap, key=lambda it: (it[0], it[1]), reverse=True)]


@dataclass
class RuleState:
    rule: Rule
    filters: list[tuple[str, Callable[[Candidate], bool]]]
    bucket_of: Callable[[Candidate], str | None]
    requested: dict[str, int]
    rank_key: Callable[[Candidate], tuple[int, ...]]
    counts: dict[str, int] = field(default_factory=dict)
    pools: dict[str, Pool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.counts = {"seen": 0} | {name: 0 for name, _ in self.filters}
        for bucket, n in self.requested.items():
            self.pools[bucket] = Pool(capacity=n + max(POOL_SLACK, 10 * n))

    def offer(self, candidate: Candidate) -> None:
        self.counts["seen"] += 1
        for name, predicate in self.filters:
            if not predicate(candidate):
                return
            self.counts[name] += 1
        bucket = self.bucket_of(candidate)
        if bucket is not None and bucket in self.pools:
            self.pools[bucket].push(self.rank_key(candidate), candidate)


class Deduplicator:
    """Suite-wide positional deduplication on REF spans."""

    def __init__(self, min_distance_bp: int):
        self.min_distance_bp = min_distance_bp
        self._accepted: dict[str, list[tuple[int, int, str, str]]] = {}

    def conflict(self, chrom: str, pos: int, ref: str, alt: str) -> str | None:
        start, end = ref_span0(pos, ref)
        for a_start, a_end, a_id, a_rule in self._accepted.get(chrom, []):
            if a_id.split(":", 1)[1] == f"{chrom}:{pos}:{ref}:{alt}":
                return f"same variant already selected as {a_id} by rule {a_rule}"
            gap = max(start - a_end, a_start - end)
            if gap < 0:
                return f"overlaps {a_id} (rule {a_rule})"
            if gap < self.min_distance_bp:
                return (
                    f"{gap} bp from {a_id} (rule {a_rule}); min_distance_bp={self.min_distance_bp}"
                )
        return None

    def accept(self, variant: TestVariant) -> None:
        start, end = ref_span0(variant.pos, variant.ref)
        self._accepted.setdefault(variant.chrom, []).append(
            (start, end, variant.variant_id, variant.rule_id)
        )


class _Context:
    def __init__(self, suite: LoadedSuite, reference: Reference):
        self.suite = suite
        self.cfg = suite.config
        self.reference = reference
        self.seed = suite.config.suite.seed
        self._beds: dict[Path, IntervalSet] = {}

    def bed(self, path: Path) -> IntervalSet:
        if path not in self._beds:
            self._beds[path] = IntervalSet(path)
        return self._beds[path]

    def stratum(self, name: str) -> IntervalSet:
        return self.bed(self.suite.resolve(self.cfg.sources.strata[name]))

    def strata_of(self, chrom: str, pos: int, ref: str) -> list[str]:
        start, end = ref_span0(pos, ref)
        return [
            name
            for name in sorted(self.cfg.sources.strata)
            if self.stratum(name).overlaps(chrom, start, end)
        ]

    def ladder_bin(self, c: Candidate) -> str | None:
        if c.variant_type not in INDEL_TYPES:
            return None
        for lo, hi in self.cfg.constraints.indel_ladder:
            if lo <= c.event_length <= hi:
                return f"{lo}-{hi}bp"
        return None


def _in_bed(ctx: _Context, bed: IntervalSet, *, contain: bool) -> Callable[[Candidate], bool]:
    def predicate(c: Candidate) -> bool:
        start, end = ref_span0(c.pos, c.ref)
        return bed.contains(c.chrom, start, end) if contain else bed.overlaps(c.chrom, start, end)

    return predicate


def _common_filters(
    ctx: _Context, rule: ClinGenRule | ClinVarRule | GiabRule | ExplicitRule
) -> list[tuple[str, Callable[[Candidate], bool]]]:
    constraints = ctx.cfg.constraints
    allowed = set(constraints.allowed_variant_types)
    if not isinstance(rule, ExplicitRule) and rule.variant_types is not None:
        allowed &= set(rule.variant_types)
    filters: list[tuple[str, Callable[[Candidate], bool]]] = [
        ("variant_type", lambda c: c.variant_type in allowed)
    ]
    max_len = constraints.max_indel_length
    if max_len is not None:
        filters.append(
            (
                "max_indel_length",
                lambda c: c.variant_type not in INDEL_TYPES or c.event_length <= max_len,
            )
        )
    return filters


def _region_filters(
    ctx: _Context, rule: ClinGenRule | ClinVarRule | GiabRule
) -> list[tuple[str, Callable[[Candidate], bool]]]:
    filters: list[tuple[str, Callable[[Candidate], bool]]] = []
    if rule.region_bed is not None:
        bed = ctx.bed(ctx.suite.resolve(rule.region_bed))
        filters.append(("region_bed", _in_bed(ctx, bed, contain=False)))
    if rule.strata:
        include = [_in_bed(ctx, ctx.stratum(n), contain=False) for n in rule.strata]
        filters.append(("strata", lambda c: any(p(c) for p in include)))
    if rule.exclude_strata:
        exclude = [_in_bed(ctx, ctx.stratum(n), contain=False) for n in rule.exclude_strata]
        filters.append(("exclude_strata", lambda c: not any(p(c) for p in exclude)))
    if rule.ladder_per_bin is not None:
        filters.append(("ladder_bin", lambda c: ctx.ladder_bin(c) is not None))
    return filters


def _ladder_requested(ctx: _Context, per_bin: int) -> dict[str, int]:
    return {f"{lo}-{hi}bp": per_bin for lo, hi in ctx.cfg.constraints.indel_ladder}


def _clinvar_state(ctx: _Context, rule: ClinVarRule) -> RuleState:
    allowed_sig = {s.lower() for s in rule.clnsig_allow}
    filters: list[tuple[str, Callable[[Candidate], bool]]] = []
    if rule.variation_ids is not None:
        ids = set(rule.variation_ids)
        filters.append(("variation_ids", lambda c: c.source.record_id in ids))
    if rule.rsids is not None:
        rsids = set(rule.rsids)
        filters.append(("rsids", lambda c: c.source.rsid in rsids))
    filters += _common_filters(ctx, rule)

    def clinvar(c: Candidate) -> ClinVarAnnotation:
        assert c.clinvar is not None
        return c.clinvar

    filters.append(("clnsig_allow", lambda c: passes_allow_list(clinvar(c).tokens, allowed_sig)))
    filters.append(("min_stars", lambda c: clinvar(c).stars >= rule.min_stars))
    if rule.clnvc is not None:
        classes = {x.lower() for x in rule.clnvc}
        filters.append(("clnvc", lambda c: (clinvar(c).clnvc or "").lower() in classes))
    if rule.genes is not None:
        genes = set(rule.genes)
        filters.append(("genes", lambda c: any(g in genes for g in clinvar(c).genes)))
    if rule.conditions is not None:
        wanted = set(rule.conditions)
        filters.append(("conditions", lambda c: bool(clinvar(c).condition_ids & wanted)))
    filters += _region_filters(ctx, rule)

    bucket_of: Callable[[Candidate], str | None]
    if rule.per_gene is not None:
        assert rule.genes is not None
        rule_genes = rule.genes
        requested = {g: rule.per_gene for g in rule_genes}

        def bucket_of(c: Candidate) -> str | None:
            return next((g for g in clinvar(c).genes if g in rule_genes), None)

    elif rule.ladder_per_bin is not None:
        requested = _ladder_requested(ctx, rule.ladder_per_bin)
        bucket_of = ctx.ladder_bin
    elif rule.variation_ids is not None:
        requested = {f"VariationID {i}": 1 for i in rule.variation_ids}

        def bucket_of(c: Candidate) -> str | None:
            return f"VariationID {c.source.record_id}"

    elif rule.rsids is not None:
        requested = {rs: 1 for rs in rule.rsids}

        def bucket_of(c: Candidate) -> str | None:
            return c.source.rsid

    else:
        assert rule.count is not None
        requested = {"total": rule.count}

        def bucket_of(c: Candidate) -> str | None:
            return "total"

    def rank_key(c: Candidate) -> tuple[int, ...]:
        return (-c.stars, derive_int(ctx.seed, "rank", c.variant_id(rule.category)))

    return RuleState(rule, filters, bucket_of, requested, rank_key)


def _erepo_state(ctx: _Context, rule: ClinGenRule) -> RuleState:
    def erepo(c: Candidate) -> ErepoAllele:
        assert c.erepo is not None
        return c.erepo

    filters: list[tuple[str, Callable[[Candidate], bool]]] = []
    if rule.erepo_ids is not None:
        ids = set(rule.erepo_ids)
        filters.append(("erepo_ids", lambda c: bool(ids.intersection(erepo(c).uuids))))
    filters += _common_filters(ctx, rule)
    allowed = set(rule.classification_allow)
    # strict: every merged eRepo record must carry an allowed classification
    filters.append(("classification_allow", lambda c: erepo(c).classifications <= allowed))
    if rule.genes is not None:
        genes = set(rule.genes)
        filters.append(("genes", lambda c: any(g in genes for g in erepo(c).genes)))
    if rule.vceps is not None:
        vceps = set(rule.vceps)
        filters.append(("vceps", lambda c: any(v in vceps for v in erepo(c).vceps)))
    if rule.conditions is not None:
        wanted = set(rule.conditions)
        filters.append(("conditions", lambda c: bool(erepo(c).condition_ids & wanted)))
    if rule.held_out == "exclude":
        filters.append(("held_out_exclude", lambda c: not erepo(c).held_out))
    elif rule.held_out == "only":
        filters.append(("held_out_only", lambda c: erepo(c).held_out))
    filters += _region_filters(ctx, rule)

    bucket_of: Callable[[Candidate], str | None]
    if rule.per_gene is not None:
        assert rule.genes is not None
        requested = {g: rule.per_gene for g in rule.genes}

        def bucket_of(c: Candidate) -> str | None:
            return c.gene

    elif rule.ladder_per_bin is not None:
        requested = _ladder_requested(ctx, rule.ladder_per_bin)
        bucket_of = ctx.ladder_bin
    elif rule.erepo_ids is not None:
        requested = {f"eRepo {i}": 1 for i in rule.erepo_ids}
        wanted_ids = rule.erepo_ids

        def bucket_of(c: Candidate) -> str | None:
            return next((f"eRepo {i}" for i in wanted_ids if i in erepo(c).uuids), None)

    else:
        assert rule.count is not None
        requested = {"total": rule.count}

        def bucket_of(c: Candidate) -> str | None:
            return "total"

    def rank_key(c: Candidate) -> tuple[int, ...]:
        return (derive_int(ctx.seed, "rank", c.variant_id(rule.category)),)

    return RuleState(rule, filters, bucket_of, requested, rank_key)


def _giab_state(ctx: _Context, rule: GiabRule) -> RuleState:
    giab = ctx.cfg.sources.giab
    assert giab is not None
    confident = ctx.bed(ctx.suite.resolve(giab.confident_bed))
    filters = [("confident_bed", _in_bed(ctx, confident, contain=True))]
    filters += _common_filters(ctx, rule)
    filters += _region_filters(ctx, rule)
    bucket_of: Callable[[Candidate], str | None]
    if rule.counts is not None:
        requested = {t.value: rule.counts[t] for t in VariantType if t in rule.counts}

        def bucket_of(c: Candidate) -> str | None:
            return c.variant_type.value

    elif rule.ladder_per_bin is not None:
        requested = _ladder_requested(ctx, rule.ladder_per_bin)
        bucket_of = ctx.ladder_bin
    else:
        assert rule.count is not None
        requested = {"total": rule.count}

        def bucket_of(c: Candidate) -> str | None:
            return "total"

    def rank_key(c: Candidate) -> tuple[int, ...]:
        return (derive_int(ctx.seed, "rank", c.variant_id(rule.category)),)

    return RuleState(rule, filters, bucket_of, requested, rank_key)


def _explicit_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}>{alt}"


def _explicit_state(ctx: _Context, rule: ExplicitRule) -> RuleState:
    filters = [("contig_in_reference", lambda c: c.chrom in ctx.reference.lengths)]
    filters += _common_filters(ctx, rule)
    order = {
        _explicit_key(v.chrom, v.pos, v.ref.upper(), v.alt.upper()): i
        for i, v in enumerate(rule.variants)
    }
    requested = dict.fromkeys(order, 1)

    def bucket_of(c: Candidate) -> str | None:
        return _explicit_key(c.chrom, c.pos, c.ref, c.alt)

    def rank_key(c: Candidate) -> tuple[int, ...]:
        return (order[_explicit_key(c.chrom, c.pos, c.ref, c.alt)],)

    return RuleState(rule, filters, bucket_of, requested, rank_key)


def _make_candidate(
    ctx: _Context,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    gene: str | None,
    source: Source,
    clinvar: ClinVarAnnotation | None,
    note: str = "",
) -> Candidate:
    return Candidate(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        variant_type=classify(ref, alt),
        event_length=event_length(ref, alt),
        gene=gene,
        source=source,
        clinvar=clinvar,
        note=note,
    )


def _stream_clinvar(ctx: _Context, states: list[RuleState], warnings: list[str]) -> None:
    src = ctx.cfg.sources.clinvar
    assert src is not None
    path = ctx.suite.resolve(src.vcf)
    skipped = {"non_simple_allele": 0, "contig_not_in_reference": 0}
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            alts = rec.alts or ()
            if len(alts) != 1 or not is_simple_allele(rec.ref) or not is_simple_allele(alts[0]):
                skipped["non_simple_allele"] += 1
                continue
            chrom = ctx.reference.resolve_contig(rec.chrom)
            if chrom is None:
                skipped["contig_not_in_reference"] += 1
                continue
            try:
                ann = parse_record(rec)
            except ValueError as exc:
                raise SelectionError(f"ClinVar {rec.chrom}:{rec.pos} (ID {rec.id}): {exc}") from exc
            ref, alt = str(rec.ref).upper(), alts[0].upper()
            source = Source(
                database="ClinVar", release=src.release, record_id=ann.variation_id, rsid=ann.rsid
            )
            base = _make_candidate(
                ctx,
                chrom,
                rec.pos,
                ref,
                alt,
                gene=ann.genes[0] if ann.genes else None,
                source=source,
                clinvar=ann,
            )
            for state in states:
                rule = state.rule
                assert isinstance(rule, ClinVarRule)
                candidate = base
                if rule.genes:
                    # the gene reported for a multi-gene record is the first one the rule targets
                    gene = next((g for g in ann.genes if g in rule.genes), base.gene)
                    if gene != base.gene:
                        candidate = replace(base, gene=gene)
                state.offer(candidate)
    for reason, n in skipped.items():
        if n:
            warnings.append(f"ClinVar: {n} records skipped ({reason})")


def _stream_erepo(ctx: _Context, states: list[RuleState], warnings: list[str]) -> None:
    src = ctx.cfg.sources.clingen_erepo
    assert src is not None
    held_out_ref: HeldOutReference | None = None
    try:
        if src.held_out_reference is not None:
            ref_cfg = src.held_out_reference
            held_out_ref = load_held_out_reference(
                ctx.suite.resolve(ref_cfg.file), ref_cfg.release, ctx.reference
            )
        else:
            warnings.append(
                f"{EREPO_DATABASE}: no held_out_reference configured; held_out is false for all "
                "eRepo variants"
            )
        loaded = load_erepo(
            ctx.suite.resolve(src.file),
            release=src.release,
            reference=ctx.reference,
            held_out_reference=held_out_ref,
        )
    except ErepoError as exc:
        raise SelectionError(str(exc)) from exc
    for reason, n in loaded.skipped.items():
        unit = "alleles" if reason == "discordant_classifications" else "records"
        warnings.append(f"{EREPO_DATABASE}: {n} {unit} skipped ({reason})")
    for allele in loaded.alleles:
        p = allele.primary
        source = Source(database=EREPO_DATABASE, release=src.release, record_id=p.uuid)
        base = _make_candidate(
            ctx,
            allele.chrom,
            allele.pos,
            allele.ref,
            allele.alt,
            gene=allele.genes[0],
            source=source,
            clinvar=None,
        )
        base = replace(base, erepo=allele)
        for state in states:
            rule = state.rule
            assert isinstance(rule, ClinGenRule)
            candidate = base
            if rule.genes:
                gene = next((g for g in rule.genes if g in allele.genes), base.gene)
                if gene != base.gene:
                    candidate = replace(base, gene=gene)
            state.offer(candidate)


def _stream_giab(ctx: _Context, states: list[RuleState], warnings: list[str]) -> None:
    giab = ctx.cfg.sources.giab
    assert giab is not None
    path = ctx.suite.resolve(giab.truth_vcf)
    skipped = {"multiallelic_or_non_simple": 0, "contig_not_in_reference": 0}
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            alts = rec.alts or ()
            if len(alts) != 1 or not is_simple_allele(rec.ref) or not is_simple_allele(alts[0]):
                skipped["multiallelic_or_non_simple"] += 1
                continue
            chrom = ctx.reference.resolve_contig(rec.chrom)
            if chrom is None:
                skipped["contig_not_in_reference"] += 1
                continue
            ref, alt = str(rec.ref).upper(), alts[0].upper()
            record_id = rec.id if rec.id else f"{chrom}:{rec.pos}:{ref}:{alt}"
            source = Source(database=giab.name, release=giab.release, record_id=str(record_id))
            candidate = _make_candidate(
                ctx, chrom, rec.pos, ref, alt, gene=None, source=source, clinvar=None
            )
            for state in states:
                state.offer(candidate)
    for reason, n in skipped.items():
        if n:
            warnings.append(f"{giab.name}: {n} records skipped ({reason})")


def _feed_explicit(ctx: _Context, state: RuleState) -> None:
    rule = state.rule
    assert isinstance(rule, ExplicitRule)
    for i, v in enumerate(rule.variants):
        ref, alt = v.ref.upper(), v.alt.upper()
        if ref == alt:
            raise SelectionError(f"rule {rule.id}: variant {v.chrom}:{v.pos} has REF == ALT")
        source = Source(
            database="suite", release=ctx.cfg.suite.version, record_id=f"{rule.id}#{i + 1}"
        )
        state.offer(
            _make_candidate(
                ctx,
                v.chrom,
                v.pos,
                ref,
                alt,
                gene=v.gene,
                source=source,
                clinvar=None,
                note=v.rationale,
            )
        )


def _to_test_variant(
    ctx: _Context, rule: Rule, c: Candidate, bucket: str, rank: int, pool_size: int
) -> TestVariant:
    clinvar_info: ClinVarInfo | None = None
    if c.clinvar is not None:
        clinvar_info = ClinVarInfo(
            clnsig=list(c.clinvar.tokens),
            review_status=c.clinvar.review_status,
            stars=c.clinvar.stars,
            conditions=[cond for cond, _ in c.clinvar.conditions],
        )
    clingen_info: ClinGenInfo | None = None
    held_out = False
    if c.erepo is not None:
        a, p = c.erepo, c.erepo.primary
        held_out = a.held_out
        clingen_info = ClinGenInfo(
            vcep="; ".join(a.vceps),
            erepo_id=p.uuid,
            other_erepo_ids=list(a.uuids[1:]),
            caid=p.caid,
            clinvar_variation_id=p.clinvar_variation_id,
            erepo_release_first_seen=a.release_first_seen,
            classification=a.classification,
            met_criteria=list(p.met_criteria),
            not_met_criteria=list(p.not_met_criteria),
            guideline=p.guideline,
            conditions=a.conditions,
        )
    if isinstance(rule, ClinGenRule):
        assert c.erepo is not None
        reason = (
            f"{rule.id}: {EREPO_DATABASE} {c.erepo.classification.value} "
            f"({'; '.join(c.erepo.vceps)}){', held out' if held_out else ''}; "
            f"bucket {bucket}; rank {rank}/{pool_size}"
        )
    elif isinstance(rule, ClinVarRule):
        assert c.clinvar is not None
        reason = (
            f"{rule.id}: ClinVar {'|'.join(c.clinvar.tokens)}, stars={c.stars} "
            f"({c.clinvar.review_status}); bucket {bucket}; rank {rank}/{pool_size}"
        )
    elif isinstance(rule, GiabRule):
        reason = (
            f"{rule.id}: {c.source.database} {c.source.release} confident-region "
            f"{c.variant_type.value}; bucket {bucket}; rank {rank}/{pool_size}"
        )
    else:
        reason = f"{rule.id}: explicit synthetic variant; {c.note}"
    return TestVariant(
        variant_id=c.variant_id(rule.category),
        rule_id=rule.id,
        chrom=c.chrom,
        pos=c.pos,
        ref=c.ref,
        alt=c.alt,
        variant_type=c.variant_type,
        indel_length=indel_length(c.ref, c.alt),
        category=rule.category,
        gene=c.gene,
        source=c.source,
        clinvar=clinvar_info,
        clingen=clingen_info,
        held_out=held_out,
        strata=ctx.strata_of(c.chrom, c.pos, c.ref),
        selection_reason=reason,
    )


def _choose(
    ctx: _Context, state: RuleState, dedup: Deduplicator
) -> tuple[list[TestVariant], RuleSummary]:
    rule = state.rule
    chosen: list[TestVariant] = []
    rejections: list[Rejection] = []
    selected: dict[str, int] = {}
    missing: list[str] = []
    for bucket, n_requested in state.requested.items():
        pool = state.pools[bucket]
        ranked = pool.ranked()
        taken = 0
        rejected_here = 0
        for rank, cand in enumerate(ranked, start=1):
            if taken >= n_requested:
                break
            vid = cand.variant_id(rule.category)
            reason = dedup.conflict(cand.chrom, cand.pos, cand.ref, cand.alt)
            if reason is not None:
                rejections.append(Rejection(variant_id=vid, reason=reason))
                log.info("rule %s: rejected %s: %s", rule.id, vid, reason)
                rejected_here += 1
                continue
            tv = _to_test_variant(ctx, rule, cand, bucket, rank, len(ranked))
            dedup.accept(tv)
            chosen.append(tv)
            taken += 1
        selected[bucket] = taken
        if taken < n_requested:
            if pool.truncated:
                raise SelectionError(
                    f"rule {rule.id}, bucket {bucket}: candidate pool of {pool.capacity} "
                    "exhausted by deduplication; increase POOL_SLACK or relax min_distance_bp"
                )
            missing.append(
                f"{bucket}: requested {n_requested}, selected {taken} "
                f"({pool.seen} candidates passed all filters, {rejected_here} rejected by "
                "deduplication)"
            )
    filter_counts = list(state.counts.items())
    summary = RuleSummary(
        rule_id=rule.id,
        required=rule.required,
        samples=list(rule.samples),
        filter_counts=filter_counts,
        requested=dict(state.requested),
        selected=selected,
        rejections=rejections,
        missing=missing,
        satisfied=not missing,
    )
    return chosen, summary


def select_variants(suite: LoadedSuite, reference: Reference) -> SelectionResult:
    """Run step 1 for a suite. Never raises for unsatisfied rules; see ``selection_failures``."""
    t0 = time.perf_counter()
    progress.info("selecting variants: %d rule(s)", len(suite.config.selection))
    ctx = _Context(suite, reference)
    warnings: list[str] = []
    states: list[RuleState] = []
    for rule in suite.config.selection:
        if isinstance(rule, ClinGenRule):
            states.append(_erepo_state(ctx, rule))
        elif isinstance(rule, ClinVarRule):
            states.append(_clinvar_state(ctx, rule))
        elif isinstance(rule, GiabRule):
            states.append(_giab_state(ctx, rule))
        else:
            states.append(_explicit_state(ctx, rule))
            if not rule.variants:
                warnings.append(f"rule {rule.id}: no variants listed")

    erepo_states = [s for s in states if isinstance(s.rule, ClinGenRule)]
    clinvar_states = [s for s in states if isinstance(s.rule, ClinVarRule)]
    giab_states = [s for s in states if isinstance(s.rule, GiabRule)]
    if erepo_states:
        _stream_erepo(ctx, erepo_states, warnings)
    if clinvar_states:
        _stream_clinvar(ctx, clinvar_states, warnings)
    if giab_states:
        _stream_giab(ctx, giab_states, warnings)
    for state in states:
        if isinstance(state.rule, ExplicitRule):
            _feed_explicit(ctx, state)

    dedup = Deduplicator(suite.config.constraints.min_distance_bp)
    variants: list[TestVariant] = []
    summaries: list[RuleSummary] = []
    for state in states:
        chosen, summary = _choose(ctx, state, dedup)
        variants.extend(chosen)
        summaries.append(summary)

    variants.sort(key=lambda v: (*reference.sort_key(v.chrom, v.pos), v.ref, v.alt, v.variant_id))
    progress.info("selected %d variants in %s", len(variants), elapsed(t0))
    return SelectionResult(
        schema_version=SCHEMA_VERSION,
        suite_name=suite.config.suite.name,
        seed=suite.config.suite.seed,
        variants=variants,
        rules=summaries,
        warnings=warnings,
    )


def variants_per_sample(result: SelectionResult) -> dict[str, list[TestVariant]]:
    samples_of = {r.rule_id: r.samples for r in result.rules}
    out: dict[str, list[TestVariant]] = {}
    for v in result.variants:
        for sample in samples_of[v.rule_id]:
            out.setdefault(sample, []).append(v)
    return out


def selection_failures(result: SelectionResult, sample_ids: Iterable[str]) -> list[str]:
    """Reasons for a non-zero exit: unsatisfied required rules and empty samples."""
    failures = [
        f"required rule {r.rule_id} is unsatisfied: {'; '.join(r.missing)}"
        for r in result.rules
        if r.required and not r.satisfied
    ]
    per_sample = variants_per_sample(result)
    failures += [
        f"sample {s} has no selected variants" for s in sample_ids if not per_sample.get(s)
    ]
    return failures
