# Test fixtures

Tiny, synthetic, hand-specified inputs. Regenerate with `python tests/fixtures/make_fixtures.py`
(fixed seed; outputs are committed). Nothing here is real genomic data; gene symbols
(`GENEA`, `GENEB`, `GENEC`, `GENEX`) are fictional.

| File | Content |
|---|---|
| `ref.fa`, `ref.fa.fai` | random sequence: chr1 10 kb, chr2 10 kb, chrX 8 kb, chrY 2 kb |
| `clinvar_mini.vcf` | 32 ClinVar-like records (contigs `1`, `2`, `X`, `MT`) |
| `giab_mini.vcf`, `giab_confident.bed` | 11 GIAB-like records (sample `HG002`) + confident regions |
| `strata_cmrg.bed` | stratification BED (covers chr1:1101-1300, chr2:5701-8100, chr2:8901-9100) |
| `suite.yaml` | fixture suite using only the files above (ClinVar + GIAB, no background) |
| `erepo_mini.tsv` | 27 ClinGen-eRepo-like rows (same columns as the tab-delimited download) |
| `erepo_v2.2.0_mini.tsv` | held-out reference: Uuid/CAID of 9 rows "present in v2.2.0" |
| `erepo_badref.tsv` | one row whose HGVS states REF `A` at chr2:9800 (FASTA: `C`) |
| `hg002_mini.vcf.gz(.tbi)`, `hg002_mini.bed` | HG002-like background (sample `HG002`, chr1/chr2) + confident BED (chr1/chr2:1-9000) |
| `suite_erepo.yaml` | eRepo P/LP vs B/LB on GENEA/GENEB/GENEX, synthetic LoF, HG002 background, two male samples; region-mode simulation (`flank_bp: 600`) |
| `sim_profile.yaml` | simulation profile for the fixtures: 2 x 100 bp, fragments 300 ± 30, 60x |

## ClinVar records exercised

| VariationID | locus | purpose |
|---|---|---|
| 100001 | 1:200 SNV | `practice_guideline` (4 stars), wins ranking; has rsID, MONDO/OMIM condition |
| 100002 | 1:450 SNV | 3 stars, 250 bp from 100001 -> rejected by `min_distance_bp` |
| 100003 / 100004 | 1:1200 del3 / 1:1201 SNV | overlapping records (2 stars each), in stratum `cmrg` |
| 100005 | 1:2000 | mixed `Pathogenic\|Uncertain_significance` -> rejected by the P/LP allow-list |
| 100006 | 1:2500 | 1 star (below `min_stars: 2`) |
| 100007, 100008 | 1:1800, 1:2800 | GENEA benign / likely benign (negative controls) |
| 100009 | 1:2300 | conflicting classifications |
| 100010 | 1:1500 | VUS |
| 100011 | 1:3400 | multi-gene `GENEINFO=GENEB:9002\|GENEZ:9099` |
| 100012, 100015 | 1:4000 ins4, 1:3600 dup (0 stars) | GENEB insertions |
| 100013, 100014 | 1:4700, 1:5500 del2 | GENEB benign |
| 100016 | 1:5000 | MNV (not an allowed variant type) |
| 100017 | 1:9950 | ALT `.` (skipped) |
| 100018 | MT:100 | contig not in the reference (skipped) |
| 100021-100027 | chr2 | indel ladder: del2, ins4, del7, ins9, del12, ins15, del25 (1 star) |
| 100028 | 2:4300 del60 | exceeds `max_indel_length: 50` |
| 100029 | 2:4500 | delins GAG>CA |
| 100035 | X:500 | inside the fixture PAR (chrX:1-1000) |
| 100031-100034 | X:2000-5000 | GENEX P, LP (del1), B, LB -> hemizygous in the male sample |

## Forced reference bases (used by tests)

- `chr2:5100` = `G`: explicit synthetic SNV G>T in `suite.yaml`.
- `chr2:9499` = `C`, `chr2:9500-9505` = `AAAAAA`, `chr2:9506` = `G`: an explicit deletion
  `chr2:9503 AA>A` is not left-aligned (bcftools norm moves it to 9499).
- `chr2:9800` = `C`: an explicit variant with REF `A` there triggers the REF-mismatch error.

## ClinGen eRepo rows exercised (`erepo_mini.tsv`)

The eRepo and background files are written by the second part of `make_fixtures.py`, which uses its
own RNG and never changes the sequence, so the files above are unaffected. HGVS accessions:
`NC_000001.11` = chr1, `NC_000002.12` = chr2, `NC_000023.11` = chrX.

| Uuid | purpose |
|---|---|
| EREPO-0001 / 0004 | GENEA P SNVs at chr1:300 / chr1:400 (99 bp apart: deduplication), both in v2.2.0 |
| EREPO-0002, 0003 | GENEA LP del2, P dup3 (held out) |
| EREPO-0005 + 0006 | one allele (chr1:5100) curated twice, P + LP: merged; 0005 in v2.2.0 -> not held out |
| EREPO-0007 + 0008 | one allele, P + B: skipped (`discordant_classifications`) |
| EREPO-0009 | retracted |
| EREPO-0010 | GRCh37 HGVS only (`no_grch38_hgvs`) |
| EREPO-0011 | uncertain breakpoints (`unsupported_hgvs`) |
| EREPO-0012 | 25 bp deletion (> `max_indel_length: 20`) |
| EREPO-0013 | VUS (rejected by P/LP and B/LB allow-lists) |
| EREPO-0014-0017 | GENEA B/LB (SNVs, insertion) |
| EREPO-0020-0022 | GENEB P/LP: delins, `g.9505del` in the chr2:9500-9505 homopolymer (VCF chr2:9499), SNV |
| EREPO-0023-0025 | GENEB B/LB |
| EREPO-0030-0033 | GENEX (chrX) P, LP, B, LB: hemizygous in male samples, outside the background BED |

## Background records (`hg002_mini.vcf.gz`)

chr1:320 and chr1:2150 (near eRepo candidates), chr1:6000 (phased, `PS`), chr1:6500 (unphased
`0/1`), chr1:7000 (`1/1`), chr1:7600 (multi-allelic `1|2`), chr1:8000 (`0/0`, skipped),
chr2:5100 (same locus as the synthetic variant: removed in `positive`), chr2:7000 (deletion),
chr2:8000 (insertion), chr2:8500 (`./.`, skipped).
