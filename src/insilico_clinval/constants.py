"""Reference-genome constants (GRCh38 only)."""

from insilico_clinval.schemas.suite import Region

# Pseudo-autosomal regions of GRCh38, 1-based closed intervals.
# Source: Genome Reference Consortium, GRCh38.p14 assembly
# (https://www.ncbi.nlm.nih.gov/grc/human), "Pseudoautosomal regions" section.
GRCH38_PAR_REGIONS: tuple[Region, ...] = (
    Region(chrom="chrX", start=10_001, end=2_781_479),
    Region(chrom="chrX", start=155_701_383, end=156_030_895),
    Region(chrom="chrY", start=10_001, end=2_781_479),
    Region(chrom="chrY", start=56_887_903, end=57_217_415),
)

SEX_CHROMOSOMES: tuple[str, str] = ("chrX", "chrY")

# RefSeq accessions of the GRCh38 primary-assembly chromosomes -> UCSC-style names.
# Source: NCBI assembly report of GRCh38.p14 (GCF_000001405.40), "RefSeq-Accn" column
# (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
# GCF_000001405.40_GRCh38.p14_assembly_report.txt). The GRCh37 versions (e.g. NC_000001.10)
# are deliberately absent; chrM (NC_012920.1) is out of scope.
GRCH38_REFSEQ_CHROMOSOMES: dict[str, str] = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
}
