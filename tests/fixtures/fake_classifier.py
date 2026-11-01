"""Fake black-box ACMG classifier for adapter tests (not a real classifier).

    python fake_classifier.py INPUT.vcf.gz OUTPUT.tsv [--fail]

Writes one row per VCF record: planted records (INFO TV_ID) whose id starts with ``curated``
or ``synthetic_novel`` get "Pathogenic" with criteria "PVS1,PM2" when the record's TV_RULE
contains ``plp`` or is the synthetic LoF, otherwise "Benign" with "BA1"; background records
get "VUS". ``--fail`` exits with status 3 (a crashing system under test).
"""

from __future__ import annotations

import sys

import pysam


def main(argv: list[str]) -> int:
    if "--fail" in argv:
        print("fake classifier: simulated crash", file=sys.stderr)
        return 3
    vcf_path, out_path = argv[0], argv[1]
    with pysam.VariantFile(vcf_path) as vcf, open(out_path, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tCLASSIFICATION\tCRITERIA\tGENE\n")
        for r in vcf:
            rule = str(r.info.get("TV_RULE", ""))
            if "TV_ID" not in r.info:
                label, criteria = "VUS", ""
            elif "plp" in rule or rule.startswith("synthetic"):
                label, criteria = "Pathogenic", "PVS1,PM2"
            else:
                label, criteria = "Benign", "BA1"
            for alt in r.alts or ():
                out.write(f"{r.chrom}\t{r.pos}\t{r.ref}\t{alt}\t{label}\t{criteria}\t\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
