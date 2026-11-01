"""Generate the tiny synthetic test fixtures (run once; outputs are committed).

    python tests/fixtures/make_fixtures.py

Writes, next to this script:

- ``ref.fa`` + ``ref.fa.fai``: random sequence (fixed seed) for chr1, chr2, chrX, chrY;
- ``clinvar_mini.vcf``: ClinVar-like VCF (ClinVar contig names ``1``, ``2``, ``X``, ``MT``);
- ``giab_mini.vcf``, ``giab_confident.bed``: GIAB-like truth set and high-confidence regions;
- ``strata_cmrg.bed``: a stratification BED.

The records below are hand-specified (position, type, length, annotation); only the REF/ALT
bases are derived from the generated sequence. Indels are made left-aligned by adjusting the
reference base after the event, so ``bcftools norm`` leaves every record unchanged.
Fixed bases used by tests (see tests/fixtures/README.md) are forced at the end.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
SEED = 20260928
CONTIGS = {"chr1": 10_000, "chr2": 10_000, "chrX": 8_000, "chrY": 2_000}
CLINVAR_CONTIG = {"chr1": "1", "chr2": "2", "chrX": "X"}
BASES = "ACGT"

rng = random.Random(SEED)
seq: dict[str, list[str]] = {c: [rng.choice(BASES) for _ in range(n)] for c, n in CONTIGS.items()}

# Forced bases (1-based): used by explicit variants in the fixture suite and by tests.
FORCED: dict[tuple[str, int], str] = {
    ("chr2", 5100): "G",  # synthetic_novel explicit SNV G>T in the fixture suite
    ("chr2", 9499): "C",  # homopolymer chr2:9500-9505 for the non-left-aligned indel test
    **{("chr2", p): "A" for p in range(9500, 9506)},
    ("chr2", 9506): "G",
    ("chr2", 9800): "C",  # REF-mismatch test uses A at this position
}


def other_base(b: str, *avoid: str) -> str:
    choices = [x for x in BASES if x != b and x not in avoid]
    return rng.choice(choices)


@dataclass
class Spec:
    chrom: str
    pos: int
    kind: str  # snv | del | ins | delins | mnv | noalt
    length: int = 1
    info: dict[str, str] | None = None
    vid: str = "."
    gt: str = "0|1"
    multiallelic: bool = False


def prepare(spec: Spec) -> None:
    """Adjust the sequence so that the event is left-aligned (minimal representation)."""
    s = seq[spec.chrom]
    a = spec.pos - 1  # 0-based anchor / first base
    # deleting s[a+1 : a+1+L] is left-aligned iff s[a] != s[a+L]
    if spec.kind == "del" and s[a] == s[a + spec.length]:
        s[a + spec.length] = other_base(s[a])


inserted: dict[tuple[str, int], str] = {}
delins_alt: dict[tuple[str, int], str] = {}
mnv_alt: dict[tuple[str, int], str] = {}
snv_alt: dict[tuple[str, int], str] = {}


def alleles(spec: Spec) -> tuple[str, str]:
    s = seq[spec.chrom]
    a = spec.pos - 1
    key = (spec.chrom, spec.pos)
    if spec.kind == "snv":
        ref = s[a]
        alt = snv_alt.setdefault(key, other_base(ref))
        return ref, alt
    if spec.kind == "del":
        return "".join(s[a : a + spec.length + 1]), s[a]
    if spec.kind == "ins":
        if key not in inserted:
            ins = [rng.choice(BASES) for _ in range(spec.length)]
            ins[-1] = other_base(s[a])  # last inserted base != anchor -> left-aligned
            inserted[key] = "".join(ins)
        return s[a], s[a] + inserted[key]
    if spec.kind == "delins":
        ref = "".join(s[a : a + spec.length])
        if key not in delins_alt:
            alt = [rng.choice(BASES) for _ in range(spec.length - 1)]
            alt[0] = other_base(ref[0], ref[-1])
            alt[-1] = other_base(ref[-1], ref[0])
            delins_alt[key] = "".join(alt)
        return ref, delins_alt[key]
    if spec.kind == "mnv":
        ref = "".join(s[a : a + spec.length])
        alt = mnv_alt.setdefault(key, "".join(other_base(b) for b in ref))
        return ref, alt
    if spec.kind == "noalt":
        return s[a], "."
    raise ValueError(spec.kind)


P, LP, PLP = "Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"
B, LB, BLB = "Benign", "Likely_benign", "Benign/Likely_benign"
S4 = "practice_guideline"
S3 = "reviewed_by_expert_panel"
S2 = "criteria_provided,_multiple_submitters,_no_conflicts"
S1 = "criteria_provided,_single_submitter"
S1C = "criteria_provided,_conflicting_classifications"
S0 = "no_assertion_criteria_provided"


def cv(
    sig: str,
    rev: str,
    gene: str,
    vc: str,
    *,
    dis: str | None = None,
    dn: str | None = None,
    rs: str | None = None,
) -> dict[str, str]:
    info = {"CLNSIG": sig, "CLNREVSTAT": rev, "CLNVC": vc, "GENEINFO": gene}
    info["CLNDISDB"] = dis or "MedGen:CN169374"
    info["CLNDN"] = dn or "not_specified"
    if rs:
        info["RS"] = rs
    return info


GA, GB, GC, GX = "GENEA:9001", "GENEB:9002|GENEZ:9099", "GENEC:9003", "GENEX:9004"
GB1 = "GENEB:9002"
DIS_A = "MONDO:MONDO:0000101,OMIM:600101"
DIS_X = "MONDO:MONDO:0000104,MedGen:C0000104"
SNV_VC, DEL_VC, INS_VC = "single_nucleotide_variant", "Deletion", "Insertion"

CLINVAR: list[Spec] = [
    # GENEA (chr1): ranking, allow-list, overlap and min-distance cases
    Spec(
        "chr1",
        200,
        "snv",
        vid="100001",
        info=cv(P, S4, GA, SNV_VC, dis=DIS_A, dn="Fixture_disease_A", rs="1001"),
    ),
    Spec(
        "chr1",
        450,
        "snv",
        vid="100002",
        info=cv(P, S3, GA, SNV_VC, dis=DIS_A, dn="Fixture_disease_A"),
    ),
    Spec("chr1", 1200, "del", 3, vid="100003", info=cv(LP, S2, GA, DEL_VC)),
    Spec("chr1", 1201, "snv", vid="100004", info=cv(PLP, S2, GA, SNV_VC)),
    Spec("chr1", 1500, "snv", vid="100010", info=cv("Uncertain_significance", S2, GA, SNV_VC)),
    Spec("chr1", 1800, "snv", vid="100007", info=cv(B, S2, GA, SNV_VC)),
    Spec("chr1", 2000, "snv", vid="100005", info=cv(f"{P}|Uncertain_significance", S2, GA, SNV_VC)),
    Spec(
        "chr1",
        2300,
        "snv",
        vid="100009",
        info=cv("Conflicting_classifications_of_pathogenicity", S1C, GA, SNV_VC),
    ),
    Spec("chr1", 2500, "snv", vid="100006", info=cv(P, S1, GA, SNV_VC)),
    Spec("chr1", 2800, "snv", vid="100008", info=cv(LB, S3, GA, SNV_VC)),
    # GENEB (chr1): multi-gene GENEINFO, indels, 0-star duplication, MNV
    Spec("chr1", 3400, "snv", vid="100011", info=cv(P, S2, GB, SNV_VC)),
    Spec("chr1", 3600, "ins", 3, vid="100015", info=cv(P, S0, GB1, "Duplication")),
    Spec("chr1", 4000, "ins", 4, vid="100012", info=cv(LP, S2, GB1, INS_VC)),
    Spec("chr1", 4700, "snv", vid="100013", info=cv(B, S2, GB1, SNV_VC)),
    Spec("chr1", 5000, "mnv", 2, vid="100016", info=cv(LP, S2, GB1, "Indel")),
    Spec("chr1", 5500, "del", 2, vid="100014", info=cv(BLB, S2, GB1, DEL_VC)),
    Spec("chr1", 9950, "noalt", vid="100017", info=cv(P, S2, GB1, SNV_VC)),
    # GENEC (chr2): indel length ladder (1 star)
    Spec("chr2", 200, "del", 2, vid="100021", info=cv(P, S1, GC, DEL_VC)),
    Spec("chr2", 800, "ins", 4, vid="100022", info=cv(LP, S1, GC, "Duplication")),
    Spec("chr2", 1400, "del", 7, vid="100023", info=cv(P, S1, GC, DEL_VC)),
    Spec("chr2", 2000, "ins", 9, vid="100024", info=cv(LP, S1, GC, INS_VC)),
    Spec("chr2", 2600, "del", 12, vid="100025", info=cv(P, S1, GC, DEL_VC)),
    Spec("chr2", 3200, "ins", 15, vid="100026", info=cv(LP, S1, GC, INS_VC)),
    Spec("chr2", 3800, "del", 25, vid="100027", info=cv(P, S1, GC, DEL_VC)),
    Spec("chr2", 4300, "del", 60, vid="100028", info=cv(P, S1, GC, DEL_VC)),
    Spec("chr2", 4500, "delins", 3, vid="100029", info=cv(LP, S1, GC, "Indel")),
    # GENEX (chrX): hemizygous in male samples; chrX:500 lies in the fixture PAR
    Spec("chrX", 500, "snv", vid="100035", info=cv(P, S1, GX, SNV_VC)),
    Spec(
        "chrX",
        2000,
        "snv",
        vid="100031",
        info=cv(P, S2, GX, SNV_VC, dis=DIS_X, dn="Fixture_X_linked_disease"),
    ),
    Spec(
        "chrX",
        3000,
        "del",
        1,
        vid="100032",
        info=cv(LP, S2, GX, DEL_VC, dis=DIS_X, dn="Fixture_X_linked_disease"),
    ),
    Spec("chrX", 4000, "snv", vid="100033", info=cv(B, S2, GX, SNV_VC)),
    Spec("chrX", 5000, "snv", vid="100034", info=cv(LB, S2, GX, SNV_VC)),
]

GIAB: list[Spec] = [
    Spec("chr1", 6400, "snv"),
    Spec("chr1", 7000, "snv"),
    Spec("chr1", 7100, "snv", gt="1|1"),
    Spec("chr1", 7800, "del", 2),
    Spec("chr1", 8500, "snv"),  # outside the confident BED
    Spec("chr1", 9200, "snv", multiallelic=True),
    Spec("chr2", 5800, "snv"),
    Spec("chr2", 6400, "del", 5),
    Spec("chr2", 7100, "ins", 3),
    Spec("chr2", 8000, "snv", gt="1|1"),
    Spec("chr2", 9000, "snv"),  # outside the confident BED
]

CONFIDENT_BED = [("chr1", 6000, 8000), ("chr1", 9000, 9500), ("chr2", 5000, 8500)]
CMRG_BED = [("chr1", 1100, 1300), ("chr2", 5700, 8100), ("chr2", 8900, 9100)]


def write_fasta() -> None:
    path = HERE / "ref.fa"
    with path.open("w") as fh:
        for chrom in CONTIGS:
            fh.write(f">{chrom}\n")
            s = "".join(seq[chrom])
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    pysam.faidx(str(path))


def write_clinvar() -> None:
    lines = [
        "##fileformat=VCFv4.1",
        "##fileDate=2026-09-01",
        "##source=ClinVar-like synthetic fixture (tests/fixtures/make_fixtures.py)",
        "##reference=GRCh38",
        "##contig=<ID=1>",
        "##contig=<ID=2>",
        "##contig=<ID=X>",
        "##contig=<ID=MT>",
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">',
        '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">',
        '##INFO=<ID=CLNVC,Number=1,Type=String,Description="Variant type">',
        '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s)">',
        '##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Disease database ids">',
        '##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease names">',
        '##INFO=<ID=RS,Number=.,Type=String,Description="dbSNP ID">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    order = {"1": 0, "2": 1, "X": 2, "MT": 3}
    rows: list[tuple[int, int, str]] = []
    for spec in CLINVAR:
        ref, alt = alleles(spec)
        assert spec.info is not None
        info = ";".join(f"{k}={v}" for k, v in spec.info.items())
        contig = CLINVAR_CONTIG[spec.chrom]
        rows.append(
            (
                order[contig],
                spec.pos,
                f"{contig}\t{spec.pos}\t{spec.vid}\t{ref}\t{alt}\t.\t.\t{info}",
            )
        )
    # A record on a contig that the reference does not contain.
    rows.append(
        (
            3,
            100,
            "MT\t100\t100018\tA\tG\t.\t.\t"
            f"CLNSIG={P};CLNREVSTAT={S4};CLNVC={SNV_VC};GENEINFO=MT-TEST:9005",
        )
    )
    rows.sort()
    (HERE / "clinvar_mini.vcf").write_text("\n".join(lines + [r for _, _, r in rows]) + "\n")


def write_giab() -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=GIAB-like synthetic fixture (tests/fixtures/make_fixtures.py)",
        *(f"##contig=<ID={c},length={n}>" for c, n in CONTIGS.items()),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002",
    ]
    for spec in GIAB:
        ref, alt = alleles(spec)
        gt = spec.gt
        if spec.multiallelic:
            alt = f"{alt},{other_base(ref, alt)}"
            gt = "1|2"
        lines.append(f"{spec.chrom}\t{spec.pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}")
    (HERE / "giab_mini.vcf").write_text("\n".join(lines) + "\n")


def write_bed(name: str, rows: list[tuple[str, int, int]]) -> None:
    (HERE / name).write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows))


# --- ClinGen eRepo export and HG002-like background ---------------
# Uses its own RNG and never modifies the sequence, so the files above stay unchanged.

erng = random.Random(SEED + 1)
ACCESSION = {"chr1": "NC_000001.11", "chr2": "NC_000002.12", "chrX": "NC_000023.11"}
EREPO_COLUMNS = [
    "#Variation",
    "ClinVar Variation Id",
    "Allele Registry Id",
    "HGVS Expressions",
    "HGNC Gene Symbol",
    "Disease",
    "Mondo Id",
    "Mode of Inheritance",
    "Assertion",
    "Applied Evidence Codes (Met)",
    "Applied Evidence Codes (Not Met)",
    "Summary of interpretation",
    "PubMed Articles",
    "Expert Panel",
    "Guideline",
    "Approval Date",
    "Published Date",
    "Retracted",
    "Evidence Repo Link",
    "Uuid",
]
GENE_INFO = {
    "GENEA": ("Fixture disease A", "MONDO:0000101", "Fixture Gene A VCEP"),
    "GENEB": ("Fixture disease B", "MONDO:0000102", "Fixture Gene B VCEP"),
    "GENEX": ("Fixture X-linked disease", "MONDO:0000104", "Fixture Gene X VCEP"),
}


def _base(chrom: str, pos: int) -> str:
    return seq[chrom][pos - 1]


def _alt_base(ref: str, *avoid: str) -> str:
    return erng.choice([b for b in BASES if b != ref and b not in avoid])


@dataclass
class ErepoSpec:
    uuid: str
    gene: str
    chrom: str
    assertion: str
    edit: str  # "snv:<pos>" or an HGVS g. edit without the accession, e.g. "900_901del"
    caid: str
    in_reference: bool = False
    met: str = "PM2_Supporting, PP3"
    retracted: bool = False
    grch37_only: bool = False
    disease: tuple[str, str] | None = None


def _hgvs(spec: ErepoSpec) -> str:
    if spec.edit.startswith("snv:"):
        pos = int(spec.edit[4:])
        key = (spec.chrom, pos)
        alt = snv_alt.setdefault(key, _alt_base(_base(spec.chrom, pos)))
        edit = f"{pos}{_base(spec.chrom, pos)}>{alt}"
    else:
        edit = spec.edit
    acc = ACCESSION[spec.chrom]
    if spec.grch37_only:
        acc = acc.rsplit(".", 1)[0] + ".10"
    return f"NM_9{spec.caid[-5:]}.1:c.100A>G, {acc}:g.{edit}"


EREPO: list[ErepoSpec] = [
    # GENEA (chr1)
    ErepoSpec(
        "EREPO-0001",
        "GENEA",
        "chr1",
        "Pathogenic",
        "snv:300",
        "CA900001",
        True,
        "PVS1, PM2_Supporting",
    ),
    ErepoSpec("EREPO-0002", "GENEA", "chr1", "Likely Pathogenic", "900_901del", "CA900002"),
    ErepoSpec("EREPO-0003", "GENEA", "chr1", "Pathogenic", "1500_1502dup", "CA900003"),
    ErepoSpec("EREPO-0004", "GENEA", "chr1", "Pathogenic", "snv:400", "CA900004", True),
    # one allele curated twice (same tier): merged; the first record is in the reference
    ErepoSpec("EREPO-0005", "GENEA", "chr1", "Pathogenic", "snv:5100", "CA900005", True),
    ErepoSpec(
        "EREPO-0006",
        "GENEA",
        "chr1",
        "Likely Pathogenic",
        "snv:5100",
        "CA900005",
        disease=("Fixture disease A2", "MONDO:0000111"),
    ),
    # one allele with discordant tiers: skipped
    ErepoSpec("EREPO-0007", "GENEA", "chr1", "Pathogenic", "snv:5700", "CA900007"),
    ErepoSpec("EREPO-0008", "GENEA", "chr1", "Benign", "snv:5700", "CA900007"),
    ErepoSpec("EREPO-0009", "GENEA", "chr1", "Pathogenic", "snv:4500", "CA900009", retracted=True),
    ErepoSpec(
        "EREPO-0010", "GENEA", "chr1", "Likely Pathogenic", "snv:6300", "CA900010", grch37_only=True
    ),
    ErepoSpec(
        "EREPO-0011", "GENEA", "chr1", "Pathogenic", "(6890_6895)_(6910_6920)del", "CA900011"
    ),
    ErepoSpec("EREPO-0012", "GENEA", "chr1", "Pathogenic", "7500_7524del", "CA900012"),
    ErepoSpec("EREPO-0013", "GENEA", "chr1", "Uncertain Significance", "snv:3900", "CA900013"),
    ErepoSpec("EREPO-0014", "GENEA", "chr1", "Benign", "snv:2100", "CA900014", True, "BA1"),
    ErepoSpec(
        "EREPO-0015",
        "GENEA",
        "chr1",
        "Likely Benign",
        "snv:2700",
        "CA900015",
        met="BP4, BS3_Supporting",
    ),
    ErepoSpec("EREPO-0016", "GENEA", "chr1", "Benign", "3300_3301insTTG", "CA900016", True, "BS1"),
    ErepoSpec("EREPO-0017", "GENEA", "chr1", "Likely Benign", "snv:6600", "CA900017", met="BP4"),
    # GENEB (chr2); 9505del is 3'-shifted in the homopolymer chr2:9500-9505 -> VCF chr2:9499
    ErepoSpec("EREPO-0020", "GENEB", "chr2", "Pathogenic", "600_602delinsTT", "CA900020", True),
    ErepoSpec("EREPO-0021", "GENEB", "chr2", "Likely Pathogenic", "9505del", "CA900021"),
    ErepoSpec("EREPO-0022", "GENEB", "chr2", "Pathogenic", "snv:1200", "CA900022"),
    ErepoSpec("EREPO-0023", "GENEB", "chr2", "Benign", "snv:1800", "CA900023", True, "BA1"),
    ErepoSpec("EREPO-0024", "GENEB", "chr2", "Likely Benign", "2400del", "CA900024", met="BP4"),
    ErepoSpec("EREPO-0025", "GENEB", "chr2", "Benign", "snv:3000", "CA900025", met="BS1, BP4"),
    # GENEX (chrX, hemizygous in male samples)
    ErepoSpec("EREPO-0030", "GENEX", "chrX", "Pathogenic", "snv:2000", "CA900030", True),
    ErepoSpec("EREPO-0031", "GENEX", "chrX", "Likely Pathogenic", "snv:3000", "CA900031"),
    ErepoSpec("EREPO-0032", "GENEX", "chrX", "Benign", "snv:4000", "CA900032", True, "BA1"),
    ErepoSpec("EREPO-0033", "GENEX", "chrX", "Likely Benign", "snv:5000", "CA900033", met="BP4"),
]
# A record whose stated reference base (A) differs from the FASTA (chr2:9800 = C).
EREPO_BADREF = ErepoSpec("EREPO-0099", "GENEB", "chr2", "Pathogenic", "9800A>G", "CA900099")


def _erepo_row(spec: ErepoSpec) -> list[str]:
    disease, mondo, vcep = GENE_INFO[spec.gene]
    if spec.disease is not None:
        disease, mondo = spec.disease
    return [
        f"{spec.gene} fixture variant {spec.uuid}",
        "",
        spec.caid,
        _hgvs(spec),
        spec.gene,
        disease,
        mondo,
        "Autosomal dominant inheritance" if spec.chrom != "chrX" else "X-linked inheritance",
        spec.assertion,
        spec.met,
        "PS4",
        "Synthetic fixture record.",
        "",
        vcep,
        f"{vcep} specifications v1",
        "2025-01-15",
        "2025-02-01",
        "true" if spec.retracted else "false",
        f"https://erepo.example.invalid/{spec.uuid}",
        spec.uuid,
    ]


def _write_tsv(name: str, header: list[str], rows: list[list[str]]) -> None:
    text = "\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows)
    (HERE / name).write_text(text)


def write_erepo() -> None:
    _write_tsv("erepo_mini.tsv", EREPO_COLUMNS, [_erepo_row(s) for s in EREPO])
    _write_tsv("erepo_badref.tsv", EREPO_COLUMNS, [_erepo_row(EREPO_BADREF)])
    # The held-out reference (an older, smaller export with fewer columns).
    ref_rows = [[s.uuid, s.caid, s.gene, s.assertion] for s in EREPO if s.in_reference]
    _write_tsv(
        "erepo_v2.2.0_mini.tsv",
        ["Uuid", "Allele Registry Id", "HGNC Gene Symbol", "Assertion"],
        ref_rows,
    )


def _left_aligned_del(chrom: str, near: int, length: int) -> int:
    pos = near
    while _base(chrom, pos) == _base(chrom, pos + length):
        pos += 1
    return pos


def write_background() -> None:
    """HG002-like background VCF (bgzip + tabix) and its high-confidence BED."""
    rows: list[tuple[str, int, str, str, str, str]] = []  # chrom, pos, ref, alt, gt, ps

    def snv(chrom: str, pos: int, gt: str, ps: str = ".") -> None:
        ref = _base(chrom, pos)
        rows.append((chrom, pos, ref, _alt_base(ref), gt, ps))

    snv("chr1", 320, "0|1")  # near EREPO-0001 / EREPO-0004
    snv("chr1", 2150, "1|0")  # near EREPO-0014
    snv("chr1", 6000, "0|1", "6000")  # phased block
    snv("chr1", 6500, "0/1")  # unphased
    snv("chr1", 7000, "1/1")
    ref = _base("chr1", 7600)
    a1 = _alt_base(ref)
    rows.append(("chr1", 7600, ref, f"{a1},{_alt_base(ref, a1)}", "1|2", "."))  # multi-allelic
    snv("chr1", 8000, "0/0")  # no alternate allele: skipped
    snv("chr2", 5100, "0|1")  # same locus as the synthetic variant of the eRepo suite
    p = _left_aligned_del("chr2", 7000, 3)
    rows.append(("chr2", p, "".join(seq["chr2"][p - 1 : p + 3]), _base("chr2", p), "0|1", "."))
    anchor = _base("chr2", 8000)
    ins = "".join(erng.choice(BASES) for _ in range(1)) + _alt_base(anchor)
    rows.append(("chr2", 8000, anchor, anchor + ins, "1|0", "."))
    snv("chr2", 8500, "./.")  # missing genotype: skipped
    lines = [
        "##fileformat=VCFv4.2",
        "##source=GIAB HG002-like synthetic fixture (tests/fixtures/make_fixtures.py)",
        *(f"##contig=<ID={c},length={n}>" for c, n in CONTIGS.items()),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002",
    ]
    for chrom, pos, r, a, gt, ps in rows:
        lines.append(f"{chrom}\t{pos}\t.\t{r}\t{a}\t50\tPASS\t.\tGT:PS\t{gt}:{ps}")
    path = HERE / "hg002_mini.vcf"
    path.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(path), preset="vcf", force=True)  # -> hg002_mini.vcf.gz (+ .tbi)
    write_bed("hg002_mini.bed", [("chr1", 0, 9000), ("chr2", 0, 9000)])


def main() -> None:
    for spec in CLINVAR + GIAB:
        prepare(spec)
    for (chrom, pos), base in FORCED.items():
        seq[chrom][pos - 1] = base
    # Derive alleles only after all sequence adjustments.
    write_fasta()
    write_clinvar()
    write_giab()
    write_bed("giab_confident.bed", CONFIDENT_BED)
    write_bed("strata_cmrg.bed", CMRG_BED)
    write_erepo()
    write_background()


if __name__ == "__main__":
    main()
