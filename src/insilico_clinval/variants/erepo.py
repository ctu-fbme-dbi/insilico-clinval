"""ClinGen Evidence Repository (eRepo) export: parsing, GRCh38 mapping, allele merging, held_out.

Input is the tab-delimited download of https://erepo.clinicalgenome.org/evrepo/ (one row per
expert-panel classification). Columns are found by header name (constants below); a missing
required column is an error that lists the header found.

Per row:

1. rows marked ``Retracted`` are skipped;
2. ``Assertion`` is mapped to the five-tier :class:`Classification` (unknown values are errors);
3. the GRCh38 genomic HGVS (``NC_0000NN.xx:g.``, see ``constants.GRCH38_REFSEQ_CHROMOSOMES``) is
   taken from ``HGVS Expressions``, converted to VCF alleles against the reference FASTA and
   left-normalised. Rows without a GRCh38 expression, with unsupported HGVS or with conflicting
   GRCh38 expressions are skipped and counted by reason. Reference bases stated in the HGVS that
   differ from the FASTA are an error (:class:`ErepoMismatchError`), reported all at once.

Rows are then grouped by normalised allele. Several rows of one allele (e.g. curated for different
conditions) are merged if their classifications are in the same tier (P/LP, VUS, B/LB); alleles
with rows in different tiers are skipped as ``discordant_classifications``.

``held_out``: with a held-out reference (e.g. eRepo v2.2.0, the release BIAS-2015 was tuned and
benchmarked on), an allele is held out when it is absent from the reference. The reference is an
earlier eRepo export (matched by ``Uuid`` or ``Allele Registry Id``) or a GRCh38 VCF of its alleles
(matched by left-normalized allele; e.g. the BIAS-2015 benchmark VCF lifted to GRCh38). Without a
reference, ``held_out`` is false and a warning is recorded.

The export is small (about 10^4 rows), so it is read into memory; it is not a ClinVar release.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from insilico_clinval.constants import GRCH38_REFSEQ_CHROMOSOMES
from insilico_clinval.genome.alleles import is_simple_allele, left_normalize
from insilico_clinval.genome.hgvs import HgvsMismatch, parse_genomic, to_vcf
from insilico_clinval.genome.reference import Reference
from insilico_clinval.schemas.enums import Classification
from insilico_clinval.schemas.variant import Condition

log = logging.getLogger(__name__)

DATABASE = "ClinGen eRepo"

COL_UUID = "Uuid"
COL_HGVS = "HGVS Expressions"
COL_GENE = "HGNC Gene Symbol"
COL_ASSERTION = "Assertion"
COL_PANEL = "Expert Panel"
COL_CAID = "Allele Registry Id"
COL_CLINVAR = "ClinVar Variation Id"
COL_DISEASE = "Disease"
COL_MONDO = "Mondo Id"
COL_MET = "Applied Evidence Codes (Met)"
COL_NOT_MET = "Applied Evidence Codes (Not Met)"
COL_GUIDELINE = "Guideline"
COL_RETRACTED = "Retracted"
COL_APPROVAL = "Approval Date"

REQUIRED_COLUMNS = (COL_UUID, COL_HGVS, COL_GENE, COL_ASSERTION, COL_PANEL)

ASSERTION_TO_CLASS: dict[str, Classification] = {
    "pathogenic": Classification.PATHOGENIC,
    "likely pathogenic": Classification.LIKELY_PATHOGENIC,
    "uncertain significance": Classification.UNCERTAIN_SIGNIFICANCE,
    "likely benign": Classification.LIKELY_BENIGN,
    "benign": Classification.BENIGN,
}

TIER: dict[Classification, str] = {
    Classification.PATHOGENIC: "P/LP",
    Classification.LIKELY_PATHOGENIC: "P/LP",
    Classification.UNCERTAIN_SIGNIFICANCE: "VUS",
    Classification.LIKELY_BENIGN: "B/LB",
    Classification.BENIGN: "B/LB",
}

_SEVERITY = list(Classification)
_TRUE = {"true", "yes", "y", "1"}
_MITOCHONDRION = "NC_012920"


class ErepoError(ValueError):
    """The eRepo export cannot be read (missing column, unknown assertion, empty Uuid …)."""


class ErepoMismatchError(ValueError):
    """Reference bases stated in eRepo HGVS differ from the reference FASTA."""


@dataclass(frozen=True)
class ErepoRecord:
    uuid: str
    caid: str | None
    clinvar_variation_id: str | None
    gene: str
    classification: Classification
    vcep: str
    guideline: str | None
    condition: Condition | None
    met_criteria: tuple[str, ...]
    not_met_criteria: tuple[str, ...]
    hgvs: str


@dataclass(frozen=True)
class ErepoAllele:
    chrom: str
    pos: int
    ref: str
    alt: str
    records: tuple[ErepoRecord, ...]  # sorted by Uuid
    held_out: bool
    release_first_seen: str

    @property
    def primary(self) -> ErepoRecord:
        return self.records[0]

    @property
    def uuids(self) -> tuple[str, ...]:
        return tuple(r.uuid for r in self.records)

    @property
    def genes(self) -> tuple[str, ...]:
        return tuple(sorted({r.gene for r in self.records}))

    @property
    def classifications(self) -> frozenset[Classification]:
        return frozenset(r.classification for r in self.records)

    @property
    def classification(self) -> Classification:
        return min(self.classifications, key=_SEVERITY.index)

    @property
    def vceps(self) -> tuple[str, ...]:
        return tuple(sorted({r.vcep for r in self.records}))

    @property
    def conditions(self) -> list[Condition]:
        out: list[Condition] = []
        for r in self.records:
            if r.condition is not None and r.condition not in out:
                out.append(r.condition)
        return out

    @property
    def condition_ids(self) -> set[str]:
        return {c.id for c in self.conditions if c.id is not None}


@dataclass
class ErepoLoad:
    alleles: list[ErepoAllele]
    rows: int
    skipped: dict[str, int] = field(default_factory=dict)
    held_out_reference: str | None = None


def _read_rows(path: Path, required: tuple[str, ...]) -> Iterator[tuple[int, dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = list(reader.fieldnames or [])
        missing = [c for c in required if c not in header]
        if missing:
            raise ErepoError(
                f"{path}: eRepo export lacks columns {missing}; found header {header}. "
                "Update the column constants in insilico_clinval/variants/erepo.py if the "
                "download format changed."
            )
        for line_no, row in enumerate(reader, start=2):
            yield line_no, {k: (v or "").strip() for k, v in row.items() if k is not None}


# An indel counts as present in a VCF reference when the reference has any indel record starting
# within this distance of the indel's span. The BIAS-2015 benchmark VCF mangles indel alleles
# (multi-base deletions written as one base, inserted sequences garbled; its erepo_to_vcf.py was
# fixed in BIAS-2015 v3, issue #39), so exact matching would wrongly hold them out. The window errs
# towards "present", i.e. towards fewer held-out variants.
INDEL_WINDOW_BP = 3


@dataclass(frozen=True)
class HeldOutReference:
    """Identifiers, alleles and indel sites of the held-out reference release."""

    release: str
    ids: frozenset[str] = frozenset()
    alleles: frozenset[tuple[str, int, str, str]] = frozenset()
    indel_sites: frozenset[tuple[str, int]] = frozenset()  # raw positions of indel records

    def contains(self, allele: tuple[str, int, str, str], ids: Iterable[str | None]) -> bool:
        if allele in self.alleles or any(i is not None and i in self.ids for i in ids):
            return True
        chrom, pos, ref, alt = allele
        if len(ref) == len(alt) or not self.indel_sites:
            return False
        window = range(pos - INDEL_WINDOW_BP, pos + len(ref) + INDEL_WINDOW_BP + 1)
        return any((chrom, p) in self.indel_sites for p in window)


def reference_ids(path: Path) -> set[str]:
    """Uuids and Allele Registry ids of an earlier eRepo export (the held-out reference)."""
    with path.open(newline="", encoding="utf-8-sig") as fh:
        header = next(csv.reader(fh, delimiter="\t"), [])
    columns = tuple(c for c in (COL_UUID, COL_CAID) if c in header)
    if not columns:
        raise ErepoError(
            f"{path}: held-out reference needs a {COL_UUID!r} or {COL_CAID!r} column; "
            f"found header {header}"
        )
    ids: set[str] = set()
    for _, row in _read_rows(path, ()):
        ids.update(row[c] for c in columns if row.get(c))
    return ids


def reference_alleles(
    path: Path, reference: Reference
) -> tuple[set[tuple[str, int, str, str]], set[tuple[str, int]]]:
    """Left-normalized GRCh38 alleles of a VCF; records on contigs absent from the FASTA or with
    a REF that differs from it are an error (wrong build). An empty ALT is a deletion of REF.
    Also returns the raw (contig, position) of every indel record."""
    alleles: set[tuple[str, int, str, str]] = set()
    indel_sites: set[tuple[str, int]] = set()
    bad: list[str] = []
    with pysam.VariantFile(str(path)) as vcf:
        for rec in vcf:
            chrom = reference.resolve_contig(rec.chrom)
            ref = (rec.ref or "").upper()
            if chrom is None or reference.fetch(chrom, rec.pos - 1, rec.pos - 1 + len(ref)) != ref:
                bad.append(f"{rec.chrom}:{rec.pos} {rec.ref}")
                continue
            for alt in rec.alts or ():
                if alt.startswith("<") or alt == "*":
                    continue
                # Every record of a reference is a variant: an empty ALT (read as ".") is an
                # unanchored deletion of REF, as in the BIAS-2015 benchmark VCF; normalization
                # adds the anchor base.
                alt = "" if alt == "." else alt
                if len(alt) != len(ref):
                    indel_sites.add((chrom, rec.pos))
                pos, r, a = left_normalize(chrom, rec.pos, ref, alt, reference)
                alleles.add((chrom, pos, r, a))
    if bad:
        raise ErepoError(
            f"{path}: {len(bad)} held-out reference records do not match the reference FASTA "
            f"(wrong build?), e.g. {bad[:3]}"
        )
    return alleles, indel_sites


def approval_dates(path: Path) -> dict[str, str]:
    """``Approval Date`` (ISO, first 10 characters) of every Uuid of an eRepo export."""
    return {
        row[COL_UUID]: row[COL_APPROVAL][:10]
        for _, row in _read_rows(path, (COL_UUID, COL_APPROVAL))
        if row.get(COL_APPROVAL)
    }


def load_held_out_reference(path: Path, release: str, reference: Reference) -> HeldOutReference:
    """A VCF (``.vcf``, ``.vcf.gz``) gives alleles; any other file is read as an eRepo export."""
    if path.name.endswith((".vcf", ".vcf.gz")):
        alleles, sites = reference_alleles(path, reference)
        return HeldOutReference(release, alleles=frozenset(alleles), indel_sites=frozenset(sites))
    return HeldOutReference(release, ids=frozenset(reference_ids(path)))


def _split_codes(value: str) -> tuple[str, ...]:
    return tuple(c.strip() for c in value.replace(";", ",").split(",") if c.strip())


def _mondo(value: str) -> str | None:
    value = value.strip()
    if not value:
        return None
    return value.replace("MONDO_", "MONDO:") if value.startswith("MONDO_") else value


def _classification(value: str, where: str) -> Classification:
    key = " ".join(value.replace("_", " ").lower().split())
    if key not in ASSERTION_TO_CLASS:
        raise ErepoError(f"{where}: unknown Assertion {value!r}")
    return ASSERTION_TO_CLASS[key]


def _grch38_allele(
    hgvs_field: str, reference: Reference
) -> tuple[tuple[str, int, str, str], str] | str:
    """Normalised ``(chrom, pos, ref, alt)`` and the HGVS used, or a skip reason."""
    expressions = [e.strip() for e in hgvs_field.split(",") if e.strip()]
    grch38 = [e for e in expressions if e.split(":", 1)[0] in GRCH38_REFSEQ_CHROMOSOMES]
    if not grch38:
        if any(e.startswith(_MITOCHONDRION) for e in expressions):
            return "mitochondrial"
        return "no_grch38_hgvs"
    alleles: dict[tuple[str, int, str, str], str] = {}
    for expr in grch38:
        chrom = reference.resolve_contig(GRCH38_REFSEQ_CHROMOSOMES[expr.split(":", 1)[0]])
        if chrom is None:
            return "contig_not_in_reference"
        try:
            h = parse_genomic(expr)
            pos, ref, alt = to_vcf(h, chrom, reference, reference.lengths[chrom])
            pos, ref, alt = left_normalize(chrom, pos, ref, alt, reference)
        except HgvsMismatch:
            raise
        except ValueError:  # HgvsUnsupported, or an indel that cannot be anchored
            return "unsupported_hgvs"
        if not (is_simple_allele(ref) and is_simple_allele(alt)):
            return "unsupported_hgvs"
        alleles.setdefault((chrom, pos, ref, alt), expr)
    if len(alleles) > 1:
        return "ambiguous_grch38_hgvs"
    ((key, expr),) = alleles.items()
    return key, expr


def load_erepo(
    path: Path,
    *,
    release: str,
    reference: Reference,
    held_out_reference: HeldOutReference | None,
) -> ErepoLoad:
    """Parse the export into merged alleles sorted by reference position."""
    skipped: dict[str, int] = {}
    mismatches: list[str] = []
    groups: dict[tuple[str, int, str, str], list[ErepoRecord]] = {}
    rows = 0

    def skip(reason: str) -> None:
        skipped[reason] = skipped.get(reason, 0) + 1

    for line_no, row in _read_rows(path, REQUIRED_COLUMNS):
        rows += 1
        where = f"{path.name}:{line_no}"
        if row.get(COL_RETRACTED, "").lower() in _TRUE:
            skip("retracted")
            continue
        uuid = row[COL_UUID]
        if not uuid:
            raise ErepoError(f"{where}: empty {COL_UUID}")
        classification = _classification(row[COL_ASSERTION], where)
        try:
            mapped = _grch38_allele(row[COL_HGVS], reference)
        except HgvsMismatch as exc:
            mismatches.append(f"{where} (Uuid {uuid}): {exc}")
            continue
        if isinstance(mapped, str):
            skip(mapped)
            continue
        key, expr = mapped
        disease, mondo = row.get(COL_DISEASE, ""), _mondo(row.get(COL_MONDO, ""))
        condition = (
            Condition(name=disease or "not provided", id=mondo) if disease or mondo else None
        )
        groups.setdefault(key, []).append(
            ErepoRecord(
                uuid=uuid,
                caid=row.get(COL_CAID) or None,
                clinvar_variation_id=row.get(COL_CLINVAR) or None,
                gene=row[COL_GENE],
                classification=classification,
                vcep=row[COL_PANEL],
                guideline=row.get(COL_GUIDELINE) or None,
                condition=condition,
                met_criteria=_split_codes(row.get(COL_MET, "")),
                not_met_criteria=_split_codes(row.get(COL_NOT_MET, "")),
                hgvs=expr,
            )
        )
    if mismatches:
        raise ErepoMismatchError(
            "eRepo HGVS reference bases differ from the reference FASTA (wrong build or "
            "reference?):\n  " + "\n  ".join(mismatches)
        )

    alleles: list[ErepoAllele] = []
    for key in sorted(groups, key=lambda k: (*reference.sort_key(k[0], k[1]), k[2], k[3])):
        records = tuple(sorted(groups[key], key=lambda r: r.uuid))
        tiers = {TIER[r.classification] for r in records}
        if len(tiers) > 1:
            skip("discordant_classifications")
            log.info(
                "eRepo allele %s:%d %s>%s skipped: discordant classifications %s",
                *key,
                {r.uuid: r.classification.value for r in records},
            )
            continue
        if held_out_reference is None:
            held_out, first_seen = False, release
        else:
            ids = [i for r in records for i in (r.uuid, r.caid)]
            present = held_out_reference.contains(key, ids)
            held_out = not present
            first_seen = held_out_reference.release if present else release
        chrom, pos, ref, alt = key
        alleles.append(ErepoAllele(chrom, pos, ref, alt, records, held_out, first_seen))
    return ErepoLoad(
        alleles=alleles,
        rows=rows,
        skipped=dict(sorted(skipped.items())),
        held_out_reference=held_out_reference.release if held_out_reference else None,
    )
