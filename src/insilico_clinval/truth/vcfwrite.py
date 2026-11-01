"""Truth VCF writing with pysam: contigs from the FASTA index, FORMAT GT (+PS) only.

Planted variants and (optionally) background records are written contig by contig in FASTA-index
order and merged by position; background records carry the INFO flag ``TV_BG``.
"""

from __future__ import annotations

import heapq
import logging
import re
from collections.abc import Iterator
from pathlib import Path

import pysam

from insilico_clinval import __version__
from insilico_clinval.genome.reference import Reference
from insilico_clinval.truth.background import BackgroundMerger, BackgroundRecord
from insilico_clinval.truth.genotype import Genotype

WRITE_LOG_EVERY = 1_000_000

log = logging.getLogger(__name__)

_UNSAFE = re.compile(r"[;=,\s]")


def _info_value(value: str) -> str:
    return _UNSAFE.sub("_", value)


def truth_header(
    sample_id: str, reference: Reference, background: BackgroundMerger | None = None
) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_line(f"##source=insilico-clinval {__version__}")
    header.add_line(f"##reference={reference.fasta.name}")
    if background is not None:
        bg = background.config
        header.add_line(
            f"##clinval_background=<Name={_info_value(bg.name)},"
            f"Release={_info_value(bg.release)},Sample={_info_value(background.column)}>"
        )
    for name, length in reference.contigs:
        header.contigs.add(name, length=length)
    for key, number, vtype, desc in (
        ("TV_ID", 1, "String", "TestVariant identifier"),
        ("TV_RULE", 1, "String", "Selection rule that selected the variant"),
        (
            "TV_CATEGORY",
            1,
            "String",
            "Variant category (standard, curated, difficult, synthetic_novel)",
        ),
        ("TV_SOURCE", 1, "String", "Source database"),
        ("TV_RELEASE", 1, "String", "Source database release"),
        ("TV_RECORD", 1, "String", "Record identifier in the source database"),
        ("TV_STARS", 1, "Integer", "ClinVar review stars (ClinVar-sourced variants only)"),
        ("TV_VCEP", 1, "String", "ClinGen expert panel (eRepo-sourced variants only)"),
        (
            "TV_HELD_OUT",
            0,
            "Flag",
            "Absent from the held-out reference release (e.g. eRepo v2.2.0)",
        ),
        ("TV_BG", 0, "Flag", "Background genome record (not a planted test variant)"),
    ):
        header.info.add(key, number, vtype, desc)
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("PS", 1, "Integer", "Phase set (position of the first member)")
    header.add_sample(sample_id)
    return header


def _key(item: Genotype | BackgroundRecord) -> tuple[int, str, tuple[str, ...], int]:
    if isinstance(item, Genotype):
        return (item.variant.pos, item.variant.ref, (item.variant.alt,), 0)
    return (item.pos, item.ref, item.alts, 1)


def write_truth_vcf(
    path: Path,
    sample_id: str,
    genotypes: list[Genotype],
    reference: Reference,
    background: BackgroundMerger | None = None,
    *,
    log_every: int = WRITE_LOG_EVERY,
) -> tuple[Path, int]:
    """Write an uncompressed, coordinate-sorted, single-sample truth VCF.

    Returns the path and the number of records; the count is logged every ``log_every`` records.
    """
    header = truth_header(sample_id, reference, background)
    by_contig: dict[str, list[Genotype]] = {}
    for g in genotypes:
        by_contig.setdefault(g.variant.chrom, []).append(g)
    n = 0
    with pysam.VariantFile(str(path), "w", header=header) as out:
        for chrom, _ in reference.contigs:
            planted: list[Genotype | BackgroundRecord] = sorted(by_contig.get(chrom, []), key=_key)
            bg: Iterator[Genotype | BackgroundRecord] = (
                background.records(chrom) if background is not None else iter(())
            )
            for item in heapq.merge(planted, bg, key=_key):
                if isinstance(item, Genotype):
                    _write_planted(out, sample_id, item)
                else:
                    _write_background(out, sample_id, item)
                n += 1
                if n % log_every == 0:
                    log.info(
                        "%s: truth VCF: %s records written (at %s)", sample_id, f"{n:,}", chrom
                    )
    return path, n


def _write_planted(out: pysam.VariantFile, sample_id: str, g: Genotype) -> None:
    v = g.variant
    rec = out.new_record(contig=v.chrom, start=v.pos - 1, alleles=(v.ref, v.alt), filter="PASS")
    rec.info["TV_ID"] = _info_value(v.variant_id)
    rec.info["TV_RULE"] = _info_value(v.rule_id)
    rec.info["TV_CATEGORY"] = v.category.value
    rec.info["TV_SOURCE"] = _info_value(v.source.database)
    rec.info["TV_RELEASE"] = _info_value(v.source.release)
    rec.info["TV_RECORD"] = _info_value(v.source.record_id)
    if v.clinvar is not None:
        rec.info["TV_STARS"] = v.clinvar.stars
    if v.clingen is not None:
        rec.info["TV_VCEP"] = _info_value(v.clingen.vcep)
    if v.held_out:
        rec.info["TV_HELD_OUT"] = True
    call = rec.samples[sample_id]
    call["GT"] = g.gt
    if len(g.gt) == 2:
        call.phased = True
    if g.phase_set is not None:
        call["PS"] = g.phase_set
    out.write(rec)


def _write_background(out: pysam.VariantFile, sample_id: str, b: BackgroundRecord) -> None:
    rec = out.new_record(contig=b.chrom, start=b.pos - 1, alleles=(b.ref, *b.alts), filter="PASS")
    rec.info["TV_BG"] = True
    call = rec.samples[sample_id]
    call["GT"] = b.gt
    if len(b.gt) == 2:
        call.phased = True
    out.write(rec)
