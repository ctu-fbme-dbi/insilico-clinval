"""Toy ACMG classifier for the quickstart (NOT a real classifier).

    python toy_classifier.py KNOWLEDGE.tsv INPUT.vcf.gz OUTPUT.tsv [--no-pp5-bp6]

Looks every called allele up in a "knowledge base" of expert classifications (here the curated
rows of the run's own ``variants.selected.tsv``, i.e. the shortcut that ACMG criteria PP5/BP6 take)
and calls everything else a VUS. A variant absent from the knowledge base, such as the synthetic
novel loss-of-function variant of the fixture suite, is therefore never reported: the quickstart
report shows it as a missed finding attributed to the interpretation stage. ``--no-pp5-bp6``
ignores the knowledge base (every variant becomes a VUS): the quickstart's "candidate" pipeline,
which ``clinval diff`` flags as a regression against the baseline.
"""

from __future__ import annotations

import csv
import sys

import pysam

LABELS = {
    "pathogenic": "Pathogenic",
    "likely_pathogenic": "Likely pathogenic",
    "likely_benign": "Likely benign",
    "benign": "Benign",
}


def main(argv: list[str]) -> int:
    use_knowledge = "--no-pp5-bp6" not in argv
    knowledge_path, vcf_path, out_path = (a for a in argv if a != "--no-pp5-bp6")
    known: dict[tuple[str, int, str, str], tuple[str, str]] = {}
    with open(knowledge_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t")) if use_knowledge else []
        for row in rows:
            if row["category"] == "curated" and row["erepo_classification"] in LABELS:
                key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
                known[key] = (LABELS[row["erepo_classification"]], row["gene"])
    with pysam.VariantFile(vcf_path) as vcf, open(out_path, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tCLASSIFICATION\tCRITERIA\tGENE\n")
        for r in vcf:
            for alt in r.alts or ():
                label, gene = known.get((r.chrom, r.pos, r.ref, alt), ("VUS", ""))
                criteria = "" if label == "VUS" else ("PP5" if "athogenic" in label else "BP6")
                out.write(f"{r.chrom}\t{r.pos}\t{r.ref}\t{alt}\t{label}\t{criteria}\t{gene}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
