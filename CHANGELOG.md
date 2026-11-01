# Changelog

## 1.0.0 (2026-11-01)

First public release; the version described in the article.

- Variant selection from ClinGen eRepo and ClinVar, synthetic novel variants, truth VCF with
  expected-outcome metadata on a GIAB HG002 background.
- Read simulation with Mason (seeded, whole genome or regions) and planting check.
- Adapters for nf-core/sarek and BIAS-2015 (Nirvana or VEP), and for black-box commands.
- Comparison with the truth (hap.py/vcfeval and per-finding attribution to the failing stage),
  metrics with Wilson intervals, acceptance criteria and regression diff.
- Fidelity comparison of simulated and real HG002 data.
- Configurations of all test sets and experiments of the article; see `docs/reproducing.md`.
