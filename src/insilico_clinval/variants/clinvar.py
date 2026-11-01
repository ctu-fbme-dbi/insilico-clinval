"""Parsing of ClinVar VCF INFO fields (CLNSIG, CLNREVSTAT, GENEINFO, CLNDISDB, CLNDN, CLNVC, RS)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from insilico_clinval.schemas.enums import Classification
from insilico_clinval.schemas.variant import Condition
from insilico_clinval.variants.stars import review_stars

# ClinVar germline classification tokens (lower case) -> five-tier class used for the
# expected classification. Combined terms map to their more severe tier.
CLNSIG_TO_CLASS: dict[str, Classification] = {
    "pathogenic": Classification.PATHOGENIC,
    "likely_pathogenic": Classification.LIKELY_PATHOGENIC,
    "pathogenic/likely_pathogenic": Classification.PATHOGENIC,
    "uncertain_significance": Classification.UNCERTAIN_SIGNIFICANCE,
    "likely_benign": Classification.LIKELY_BENIGN,
    "benign": Classification.BENIGN,
    "benign/likely_benign": Classification.BENIGN,
}

_SEVERITY = list(Classification)

# Preferred condition ontologies, most preferred first.
_CONDITION_PREFIXES = ("MONDO:", "OMIM:", "Orphanet:")

# CLNDN placeholders that do not name a condition.
_PLACEHOLDER_CONDITIONS = {"not specified", "not provided", "see cases"}


def info_str(record: Any, key: str) -> str | None:
    """INFO value as one string; pysam splits Number=. strings on commas, so re-join them."""
    if key not in record.info:
        return None
    value = record.info[key]
    if isinstance(value, tuple):
        return ",".join(str(v) for v in value)
    if value is True:
        return ""
    return str(value)


def clnsig_tokens(clnsig: str) -> list[str]:
    """Split CLNSIG into lower-case tokens on ``|`` and ``,`` (leading ``_`` stripped)."""
    tokens: list[str] = []
    for part in clnsig.replace(",", "|").split("|"):
        token = part.strip().lstrip("_").strip().lower()
        if token:
            tokens.append(token)
    return tokens


def passes_allow_list(tokens: list[str], allowed: set[str]) -> bool:
    """Strict allow-list: every token must be allowed (mixed P + VUS is rejected)."""
    return bool(tokens) and all(t in allowed for t in tokens)


def expected_class(tokens: list[str]) -> Classification | None:
    classes = [CLNSIG_TO_CLASS[t] for t in tokens if t in CLNSIG_TO_CLASS]
    if not classes:
        return None
    return min(classes, key=_SEVERITY.index)


def _normalise_condition_id(raw: str) -> str | None:
    raw = raw.strip()
    if not raw or raw in (".", "not_provided", "not_specified"):
        return None
    if raw.startswith("MONDO:MONDO:"):
        return raw[len("MONDO:") :]
    if raw.startswith("Human_Phenotype_Ontology:"):
        return raw[len("Human_Phenotype_Ontology:") :]
    return raw


def parse_conditions(clndisdb: str | None, clndn: str | None) -> list[tuple[Condition, set[str]]]:
    """Pair CLNDN names with CLNDISDB id groups; the condition id prefers MONDO, OMIM, Orphanet."""
    groups = (clndisdb or "").split("|") if clndisdb else []
    names = (clndn or "").split("|") if clndn else []
    out: list[tuple[Condition, set[str]]] = []
    for i in range(max(len(groups), len(names))):
        ids_list = [
            nid
            for nid in (
                _normalise_condition_id(x)
                for x in (groups[i].split(",") if i < len(groups) else [])
            )
            if nid is not None
        ]
        name = names[i].replace("_", " ").strip() if i < len(names) else ""
        if name.lower() in _PLACEHOLDER_CONDITIONS or (not name and not ids_list):
            continue
        primary: str | None = None
        for prefix in _CONDITION_PREFIXES:
            primary = next((x for x in ids_list if x.startswith(prefix)), None)
            if primary is not None:
                break
        if primary is None and ids_list:
            primary = ids_list[0]
        out.append((Condition(name=name or "not provided", id=primary), set(ids_list)))
    return out


def parse_genes(geneinfo: str | None) -> list[str]:
    """``BRCA1:672|NBR2:10230`` -> ``["BRCA1", "NBR2"]``."""
    if not geneinfo:
        return []
    return [g.split(":")[0] for g in geneinfo.split("|") if g]


@dataclass(frozen=True)
class ClinVarAnnotation:
    variation_id: str
    tokens: list[str]
    review_status: str
    stars: int
    clnvc: str | None
    genes: list[str]
    conditions: list[tuple[Condition, set[str]]]
    rsid: str | None

    @property
    def condition_ids(self) -> set[str]:
        return {i for _, ids in self.conditions for i in ids}


def parse_record(record: Any) -> ClinVarAnnotation:
    clnsig = info_str(record, "CLNSIG") or ""
    revstat = info_str(record, "CLNREVSTAT") or ""
    rs = info_str(record, "RS")
    rsid = f"rs{rs.split(',')[0]}" if rs else None
    return ClinVarAnnotation(
        variation_id=str(record.id) if record.id else "",
        tokens=clnsig_tokens(clnsig),
        review_status=revstat,
        stars=review_stars(revstat) if revstat else 0,
        clnvc=info_str(record, "CLNVC"),
        genes=parse_genes(info_str(record, "GENEINFO")),
        conditions=parse_conditions(info_str(record, "CLNDISDB"), info_str(record, "CLNDN")),
        rsid=rsid,
    )
