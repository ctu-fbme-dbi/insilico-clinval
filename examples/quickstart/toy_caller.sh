#!/usr/bin/env bash
# Toy secondary analysis for the quickstart (NOT a clinical pipeline):
# calls variants with bcftools on Mason's true alignments, i.e. it skips read alignment.
#   toy_caller.sh FASTQ_R1 SAMPLE FASTA OUTDIR
set -euo pipefail
r1=$1 sample=$2 fasta=$3 outdir=$4
bam="$(dirname "$r1")/${sample}.true.bam"
mkdir -p "$outdir"
bcftools mpileup -Ou -a AD,DP -f "$fasta" "$bam" \
  | bcftools call -mv -Ou \
  | bcftools norm -f "$fasta" -Oz -o "$outdir/calls.vcf.gz"
bcftools index -t "$outdir/calls.vcf.gz"
