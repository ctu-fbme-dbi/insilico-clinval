"""Parsers of tertiary-analysis outputs into CanonicalFinding lists.

``reported`` follows the pipeline's reporting policy (classification in
``report_classifications``) unless the tool states it explicitly (generic ``reported`` column).
Unknown classification labels are errors: silently mapping them would hide interpretation
failures.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

import pysam

from insilico_clinval.adapters import SutFailure
from insilico_clinval.schemas.enums import Classification
from insilico_clinval.schemas.finding import CanonicalFinding
from insilico_clinval.schemas.pipeline import ParserSpec

# BIAS-2015 v3.0.0 `acmgClassification` labels (test/data/*.bias_output.tsv).
BIAS_LABELS: dict[str, Classification] = {
    "pathogenic": Classification.PATHOGENIC,
    "likely pathogenic": Classification.LIKELY_PATHOGENIC,
    "uncertain": Classification.UNCERTAIN_SIGNIFICANCE,
    "uncertain significance": Classification.UNCERTAIN_SIGNIFICANCE,
    "likely benign": Classification.LIKELY_BENIGN,
    "benign": Classification.BENIGN,
}
BIAS_COLUMNS = ("chromosome", "position", "refAllele", "altAllele", "acmgClassification")


def _none(value: str | None) -> str | None:
    value = (value or "").strip()
    return value or None


def bias_criteria(rationale: str) -> list[str]:
    """Met ACMG codes from BIAS's rationale JSON (``{group: {code: [strength, text]}}``)."""
    if not rationale.strip():
        return []
    try:
        data = json.loads(rationale)
    except json.JSONDecodeError as exc:
        raise SutFailure(f"BIAS rationale is not JSON: {rationale[:80]!r}") from exc
    met = [
        code.upper()
        for group in data.values()
        for code, value in group.items()
        if isinstance(value, list) and value and value[0]
    ]
    return sorted(met)


Allele = tuple[str, int, str, str]


def vep_allele_map(vcf: Path) -> dict[Allele, Allele]:
    """VEP's allele representation of each biallelic record in ``vcf`` → its VCF representation.

    VEP drops the leading bases REF and ALT share and shifts POS accordingly, writing an empty
    allele as ``-`` (e.g. ``chr10:87864518 CAA>C`` becomes ``chr10:87864519 AA>-``; a split
    multiallelic ``CAAAACAAAAACA>CAAAACA`` loses all seven shared bases, and then shared trailing
    bases too: ``CATTCATTTATTTATTTATTT>CATTTATTTATTT`` becomes ``+4 CATTTATT>-``; with no shared
    first base only the shared trailing bases go, ``AAAGAAGAAGAAG>GAAGAAGAAGAAGAAGAAG`` becomes
    ``A>GAAGAAG`` at POS); BIAS-2015 in VEP mode reports alleles that way. Every prefix length,
    with and without the shared suffix, is mapped, so the result does not depend on how far VEP
    trims. Records VEP leaves unchanged map to themselves; exact VCF alleles win over trimmed forms
    of other records.
    """
    exact: dict[Allele, Allele] = {}
    trimmed: dict[Allele, Allele] = {}
    with pysam.VariantFile(str(vcf)) as fh:
        for rec in fh:
            ref = rec.ref or ""
            for alt in rec.alts or ():
                key: Allele = (rec.chrom, rec.pos, ref, alt)
                exact[key] = key
                r, a = ref, alt
                while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:  # shared trailing bases only
                    r, a = r[:-1], a[:-1]
                if (r, a) != (ref, alt):
                    trimmed.setdefault((rec.chrom, rec.pos, r, a), key)
                k = 0
                while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
                    k += 1
                    if len(ref) + len(alt) <= 2:
                        continue
                    r, a = ref[k:], alt[k:]
                    trimmed.setdefault((rec.chrom, rec.pos + k, r or "-", a or "-"), key)
                    while r and a and r[-1] == a[-1]:  # then shared trailing bases
                        r, a = r[:-1], a[:-1]
                    trimmed.setdefault((rec.chrom, rec.pos + k, r or "-", a or "-"), key)
    return {**trimmed, **exact}


def parse_bias_tsv(
    path: Path,
    *,
    sample_id: str,
    tool: str,
    tool_version: str,
    report_classifications: Iterable[Classification],
    allele_map: dict[Allele, Allele] | None = None,
) -> list[CanonicalFinding]:
    """BIAS-2015 TSV → findings; ``allele_map`` (VEP mode) restores the input VCF's alleles."""
    reportable = set(report_classifications)
    findings: list[CanonicalFinding] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in BIAS_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SutFailure(f"{path.name}: BIAS output lacks columns {missing}")
        for row in reader:
            label = row["acmgClassification"].strip().lower().replace("_", " ")
            if label not in BIAS_LABELS:
                raise SutFailure(f"{path.name}: unknown BIAS classification {label!r}")
            classification = BIAS_LABELS[label]
            allele: Allele = (
                row["chromosome"],
                int(row["position"]),
                row["refAllele"].upper(),
                row["altAllele"].upper(),
            )
            if allele_map is not None:
                if allele not in allele_map:
                    raise SutFailure(f"{path.name}: {allele} matches no record of the input VCF")
                allele = allele_map[allele]
            chrom, pos, ref, alt = allele
            findings.append(
                CanonicalFinding(
                    sample_id=sample_id,
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    gene=_none(row.get("geneName")),
                    hgvs_c=_none(row.get("hgvsc")),
                    hgvs_p=_none(row.get("hgvsp")),
                    classification=classification,
                    reported=classification in reportable,
                    acmg_criteria=bias_criteria(row.get("rationale") or ""),
                    tool=tool,
                    tool_version=tool_version,
                )
            )
    return findings


def _records(path: Path, spec: ParserSpec) -> Iterator[dict[str, Any]]:
    if spec.format == "json":
        data = json.loads(path.read_text())
        if spec.records is not None:
            data = data[spec.records]
        if not isinstance(data, list):
            raise SutFailure(f"{path.name}: expected a list of records")
        yield from data
        return
    delimiter = "\t" if spec.format == "tsv" else ","
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        missing = [c for c in spec.columns.values() if c not in (reader.fieldnames or [])]
        if missing:
            raise SutFailure(f"{path.name}: missing columns {missing}")
        yield from reader


def parse_generic(
    path: Path,
    spec: ParserSpec,
    *,
    sample_id: str,
    tool: str,
    tool_version: str,
    report_classifications: Iterable[Classification],
) -> list[CanonicalFinding]:
    reportable = set(report_classifications)
    cols = spec.columns
    findings: list[CanonicalFinding] = []
    for n, record in enumerate(_records(path, spec), start=1):

        def get(key: str, record: dict[str, Any] = record) -> str | None:
            column = cols.get(key)  # type: ignore[call-overload]
            if column is None:
                return None
            value = record.get(column)
            return None if value is None else str(value).strip()

        label = get("classification")
        classification: Classification | None = None
        if label:
            if label not in spec.classification_map:
                raise SutFailure(
                    f"{path.name} record {n}: classification {label!r} not in classification_map"
                )
            classification = spec.classification_map[label]
        criteria_text = get("criteria") or ""
        criteria = sorted(
            {c.strip().upper() for c in criteria_text.split(spec.criteria_separator) if c.strip()}
        )
        reported_text = get("reported")
        if reported_text is not None:
            reported = reported_text.lower() in {v.lower() for v in spec.reported_true}
        else:
            reported = classification in reportable
        chrom, pos, ref, alt = get("chrom"), get("pos"), get("ref"), get("alt")
        if not (chrom and pos and ref and alt):
            raise SutFailure(f"{path.name} record {n}: missing chrom/pos/ref/alt")
        findings.append(
            CanonicalFinding(
                sample_id=sample_id,
                chrom=chrom,
                pos=int(pos),
                ref=ref.upper(),
                alt=alt.upper(),
                gene=get("gene") or None,
                hgvs_c=get("hgvs_c") or None,
                hgvs_p=get("hgvs_p") or None,
                classification=classification,
                reported=reported,
                acmg_criteria=criteria,
                tool=tool,
                tool_version=tool_version,
            )
        )
    return findings
