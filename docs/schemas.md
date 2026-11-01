# Data contracts

Source of truth for all data exchanged between modules. Implemented as pydantic v2 models in
`src/insilico_clinval/schemas/` (`model_config = ConfigDict(extra="forbid", frozen=True)`),
exported as JSON Schema into `docs/json_schema/` by `clinval export-schemas`.
Every persisted file carries `schema_version`.

Enumerations are lower-case strings.

## Enumerations

- `VariantCategory`: `standard` | `curated` | `difficult` | `synthetic_novel`
- `VariantType`: `snv` | `mnv` | `insertion` | `deletion` | `delins` (larger types reserved: `cnv`, `sv`, `str`)
- `Zygosity`: `het` | `hom` | `hemi`
- `Haplotype`: `1` | `2` | `both`
- `Classification` (ACMG/AMP five-tier): `pathogenic` | `likely_pathogenic` | `uncertain_significance` | `likely_benign` | `benign`
- `SampleRole`: `positive` | `negative` | `mixed`
- `Sex`: `female` (XX) | `male` (XY)
- `EvaluationLevel`: `calling` | `interpretation` | `both`
- `Phase`: `cis` | `trans` (compound heterozygotes)
- `OutcomeClass` (finding level, see below)
- `PlantingStatus`: `ok` | `missing_on_haplotype` | `no_reads`
- `SimulationMode`: `region` | `genome`

## TestVariant

One selected variant with full provenance.

| Field | Type | Notes |
|---|---|---|
| `variant_id` | str | stable, deterministic: `{category}:{chrom}:{pos}:{ref}:{alt}` |
| `rule_id` | str | selection rule that selected it |
| `chrom`, `pos`, `ref`, `alt` | str, int (1-based), str, str | GRCh38, `chr` prefix, normalised (left-aligned, biallelic) |
| `variant_type` | VariantType | |
| `indel_length` | int | `len(alt) - len(ref)`; 0 for SNV/MNV |
| `category` | VariantCategory | |
| `gene` | str \| null | HGNC symbol |
| `source` | `{database, release, record_id, rsid}` | e.g. `ClinVar`, `2026-09-01`, VariationID, `rs…`; `record_id` is always set (eRepo: `database=ClinGen eRepo`, `record_id` = Uuid of the primary record; GIAB records without an ID: `chrom:pos:ref:alt`; explicit variants: `database=suite`, `release=` suite version, `record_id={rule_id}#{n}`) |
| `clinvar` | `{clnsig: list[str], review_status: str, stars: int, conditions: list[{name, id}]}` \| null | `clnsig`: lower-case tokens; stars: see *ClinVar review stars* below |
| `strata` | list[str] | names of difficult-region BEDs the variant falls in |
| `clingen` | `{vcep, erepo_id, other_erepo_ids, caid, clinvar_variation_id, erepo_release_first_seen, classification, met_criteria, not_met_criteria, guideline, conditions: list[{name, id}]}` \| null | ClinGen eRepo allele, see *ClinGen eRepo* below |
| `held_out` | bool | true if absent from the held-out reference release (e.g. eRepo v2.2.0, on which BIAS-2015 was tuned); false for non-eRepo variants and when no reference is configured |
| `selection_reason` | str | human-readable rule that selected it |

### ClinVar review stars

Explicit mapping of `CLNREVSTAT` (`insilico_clinval/variants/stars.py`, after ClinVar's
"Review status" documentation). Unknown values are an error.

| CLNREVSTAT | stars |
|---|---:|
| `practice_guideline` | 4 |
| `reviewed_by_expert_panel` | 3 |
| `criteria_provided,_multiple_submitters,_no_conflicts` | 2 |
| `criteria_provided,_single_submitter` | 1 |
| `criteria_provided,_conflicting_classifications` (older: `…_interpretations`) | 1 |
| `no_assertion_criteria_provided`, `no_classification_provided`, `no_classification_for_the_single_variant`, `no_interpretation_for_the_single_variant`, `no_classifications_from_unflagged_records`, `no_assertion_provided` | 0 |

Conflicting classifications get **1 star**, as in ClinVar's review-status table. Such records are rejected anyway by any allow-list, because their `CLNSIG` is
`Conflicting_classifications_of_pathogenicity`.

### CLNSIG allow-list

`CLNSIG` is split on `|` and `,` into lower-case tokens (leading `_` removed). A rule's
`clnsig_allow` is **strict**: every token must be allowed, so `Pathogenic|Uncertain_significance`
and `Pathogenic|risk_factor` are rejected by a P/LP list. The expected classification is the most
severe five-tier class among the tokens (`pathogenic/likely_pathogenic` → `pathogenic`,
`benign/likely_benign` → `benign`).

### ClinGen eRepo

Source: the tab-delimited download of the ClinGen Evidence Repository
(`insilico_clinval/variants/erepo.py`). Columns are located by header name; required:
`Uuid`, `HGVS Expressions`, `HGNC Gene Symbol`, `Assertion`, `Expert Panel`; optional:
`Allele Registry Id`, `ClinVar Variation Id`, `Disease`, `Mondo Id`,
`Applied Evidence Codes (Met)`, `Applied Evidence Codes (Not Met)`, `Guideline`, `Retracted`.
A missing required column is an error that prints the header found.

- `Retracted` rows are skipped. `Assertion` maps to the five-tier class (`Pathogenic`,
  `Likely Pathogenic`, `Uncertain Significance`, `Likely Benign`, `Benign`; case and `_`
  insensitive); any other value is an error.
- Coordinates: the GRCh38 genomic HGVS (`NC_0000NN.xx:g.`, accessions in
  `constants.GRCH38_REFSEQ_CHROMOSOMES`) is converted to VCF against the reference FASTA and
  left-normalised (HGVS is 3'-shifted). Supported edits: substitution, `del`, `dup`, `ins`,
  `delins`. Rows without a GRCh38 expression (`no_grch38_hgvs`), on chrM, with unsupported
  syntax (uncertain positions, repeats, inversions …) or with conflicting GRCh38 expressions are
  skipped and counted by reason (warnings, `selection_report.md`). A reference base stated in the
  HGVS that differs from the FASTA is an error (exit 3), reported for all rows at once.
- Rows of the same normalised allele are merged when their classifications are in one tier
  (P/LP, VUS or B/LB): `erepo_id` = smallest Uuid, `classification` = most severe class, VCEPs and
  conditions are united, `met_criteria` come from the primary record. Alleles with rows in
  different tiers are skipped (`discordant_classifications`).
- `held_out`: with `sources.clingen_erepo.held_out_reference`, an allele is held out when it is
  absent from the reference. The reference is an export of the reference release (needs a `Uuid`
  or `Allele Registry Id` column; matched by Uuid or CAID) or a GRCh38 VCF of its alleles (`.vcf`,
  `.vcf.gz`; SNVs matched by left-normalized allele, indels also by any reference indel record
  starting within 3 bp of their span, because the BIAS-2015 benchmark VCF mangles indel alleles;
  an empty ALT is an unanchored deletion; a REF that differs from the FASTA is an error).
  The article uses the BIAS-2015 benchmark of eRepo V2.2.0 lifted to GRCh38
  (`bias2015_erepo_benchmark_v2.2.0` in `scripts/data_sources.yaml`). `erepo_release_first_seen` is the reference release
  for alleles present in it, else the current release. Without a reference, `held_out` is false
  for all variants, a warning is recorded, and rules with `held_out: exclude|only` are invalid.
  For an existing run, `clinval held-out` recomputes the flag without touching the selection and
  writes `held_out.tsv` (`variant_id`, `gene`, `erepo_id`, `held_out`, `first_approval`,
  `approved_after_cutoff`, `agrees_with_date`); `clinval compare` adds `held_out` to the evidence
  and to `outcomes.tsv`. The approval-date columns are a plausibility check: with
  `--approved-after` (the reference's snapshot date), an allele whose eRepo records were all
  approved after it should be held out; disagreements are listed as a manifest warning.
- `classification_allow` is strict: every merged record must carry an allowed class.

## ExpectedOutcome

What the system under test is expected to produce for one planted variant.

| Field | Type | Notes |
|---|---|---|
| `variant_id` | str | → TestVariant |
| `expected_classification` | Classification \| null | null = not assessed (e.g. pure calling test) |
| `acceptable_classifications` | list[Classification] | tolerance, e.g. P accepts LP |
| `expected_reportable` | bool | must appear as a reported finding |
| `reportable_basis` | str | why (gene/condition in scope, negative control, …) |
| `condition` | `{name, id}` \| null | MONDO/OMIM id |
| `expected_acmg_criteria` | list[str] | optional, e.g. `["PVS1","PM2"]` |
| `evaluation_level` | `calling` \| `interpretation` \| `both` | which layers are scored |
| `primary_criterion` | `reportability` \| `classification` | default `reportability` (2-tier: reported yes/no); 5-tier classification is secondary |
| `notes` | str | e.g. hemizygous conversion, compound-het phase, synthetic rationale |

Persisted as `<sample_id>.expected_outcomes.json` =
`ExpectedOutcomeSet {schema_version, sample_id, outcomes: list[ExpectedOutcome]}` (sorted by
`variant_id`). Derivation from the variant's rule `expected` block: `level: calling` ⇒
`expected_classification = null`, `acceptable_classifications = []`; otherwise the expected
classification comes from ClinVar `CLNSIG`, from the eRepo classification or from the explicit
variant's `classification` and must be one of the rule's acceptable classifications (else
error). For eRepo variants `expected_acmg_criteria` are the met criteria of the primary record,
the condition is the first eRepo condition matching the rule's `conditions` (MONDO), else the
first one, and `notes` state the eRepo id, VCEP and held-out status. `reportable_basis` defaults to
a generated sentence (e.g. `negative control: benign variant must not be reported (rule …)`).
The condition is the first ClinVar condition matching the rule's `conditions` filter, else the
first with an id (MONDO > OMIM > Orphanet preferred); CLNDN placeholders (`not specified`,
`not provided`) are skipped.

## SampleSpec (one simulated genome)

| Field | Type | Notes |
|---|---|---|
| `sample_id` | str | |
| `role` | SampleRole | positive (P/LP), negative (B/LB), mixed |
| `sex` | Sex | drives chrX/chrY ploidy (hemizygous in males); HG002 background ⇒ male |
| `background` | `{name, release, vcf, bed, records, removed, min_distance_bp, removed_log, planted_outside_bed, skipped, spanning_deletion_to_ref}` \| null | e.g. GIAB HG002 v4.2.1; `skipped` = records skipped by reason (`non_simple_allele`, `outside_bed`, `not_pass`, `missing_genotype`, `no_alt_allele`, `too_long`), `spanning_deletion_to_ref` = `*` alleles read as REF; `records` = background records written (before splitting), `removed` = records dropped within `min_distance_bp` of planted variants (the gap used is in `min_distance_bp`), listed in `removed_log` (`<sample>.background_removed.tsv`: `chrom, pos, ref, alt, gt, near_variant_id, gap_bp`); `planted_outside_bed` = planted variant ids outside the background BED (allowed, flagged) |
| `schema_version` | int | |
| `variants` | list[`{variant_id, zygosity, haplotype, phase_set}`] | `phase_set` (int \| null, = VCF `PS` = position of the first member) groups compound heterozygotes; cis = same haplotype, trans = different. Haplotype `1` ⇔ GT `1\|0`, `2` ⇔ `0\|1`, `both` ⇔ `1\|1`; hemizygous variants have haplotype `1` |
| `truth_vcf` | path | relative to the run directory; phased GT (`0\|1`, `1\|0`, `1\|1`, hemizygous `1`) |
| `expected_outcomes` | path | relative to the run directory; ExpectedOutcomeSet |

Genotype rules: zygosity = sample `overrides` > rule `zygosity` > sample `default_zygosity`.
In `male` samples chrX/chrY outside the PARs are haploid (`hemi`, GT `1`); chrY in a `female`
sample, `hemi` on a diploid locus and compound hets on hom/hemi loci are errors. Unphased
het haplotypes and the first haplotype of each phase set are drawn with `derive_choice` (see
`docs/architecture.md`, Determinism).

Background merge (`truth/background.py`): the background VCF (bgzipped + tabix) is streamed per
contig; records without an ALT allele in the background sample, with a missing genotype or with
non-simple alleles are skipped (counted in warnings); records whose REF span overlaps or lies
closer than `background.min_distance_bp` (default: `constraints.min_distance_bp`) to a planted
variant are removed and logged; unphased background
hets are phased with `derive_choice(seed, "bg_haplotype", "{sample}|{chrom}:{pos}")`; the
background `PS` is dropped and ALT alleles absent from the genotype are removed. The sample's sex
must equal the background's (`background.sex`; HG002 is male). GIAB v4.2.1 covers chr1–22 only,
so a male HG002 sample has no background truth on chrX/chrY.

Optional record filters of `BackgroundConfig` (all off by default, so the truth VCFs of
the v4.2.1 suites do not change): `restrict_to_bed` skips records whose REF span is not inside the
BED; `require_pass` skips records with a FILTER other than `PASS` or missing (`.`; GIAB v5.0q
writes `.` for passing records and `GAP1`/`GAP2`/`HET1`/… otherwise); `max_length` skips records
with a genotyped allele longer than the limit (`event_length`); `spanning_deletion: ref` reads a
`*` allele as REF on its haplotype instead of skipping the whole record (the deletion is carried
by another record). The assembly-based GIAB v5.0q background needs all four (records outside its
BED can lack one haplotype, `0|.`).

### Truth VCF

Single sample (column = `sample_id`), coordinate-sorted in FASTA-index order, bgzipped +
tabix. Header: `##contig` for every FASTA-index contig exactly once, INFO `TV_ID`, `TV_RULE`,
`TV_CATEGORY`, `TV_SOURCE`, `TV_RELEASE`, `TV_RECORD`, `TV_STARS` (ClinVar only), `TV_VCEP`
(eRepo only), `TV_HELD_OUT` (flag), `TV_BG` (flag: background record), and
`##clinval_background=<Name,Release,Sample>` when a background is merged; FORMAT `GT` and `PS`
only (no AD/DP/GQ/PL). FILTER `PASS`. Planted variants and background records are merged by
position. Written by pysam, then `bcftools norm --no-version -f REF -c e -m -any`; if
normalisation changes any planted record, the run fails (selection and variant ids are computed
on source coordinates). Background records may be split (multi-allelic) or realigned; the output
must stay coordinate-sorted.

## SuiteConfig (YAML input)

See `configs/suites/example.yaml`, `tests/fixtures/suite.yaml` and
`tests/fixtures/suite_erepo.yaml`. Unknown and duplicate keys are
errors. Relative paths are resolved against the directory of the suite file; `--out` overrides
`output.directory`.

| Key | Content |
|---|---|
| `schema_version` | `1` |
| `suite` | `name`, `version`, `seed` (int) |
| `reference` | `build: GRCh38`, `fasta` (with `.fai`), `sha256` (optional; checked if set), `par_regions` (optional list of `{chrom, start, end}`, 1-based closed; default GRCh38 PAR1/PAR2) |
| `sources` | `clingen_erepo {file, release, held_out_reference {release, file}?}`, `clinvar {vcf, release}`, `giab {truth_vcf, confident_bed, release, name=GIAB}`, `strata {name: bed}` |
| `constraints` | `min_distance_bp`, `allowed_variant_types`, `max_indel_length` (optional), `indel_ladder` (ascending disjoint `[lo, hi]` bins, bp) |
| `selection` | list of rules (below), processed in order |
| `background` | optional `{name, vcf, bed, release, sex, sample?, min_distance_bp?, restrict_to_bed?, require_pass?, max_length?, spanning_deletion?}` (filters: see Background merge; `min_distance_bp` = gap to planted variants, default `constraints.min_distance_bp`): background genome merged into every sample unless the sample sets `background: false`; `vcf` bgzipped + tabix; `sample` = VCF column (default: first) |
| `samples` | `sample_id`, `role`, `sex`, `default_zygosity` (`het`\|`hom`), `overrides [{variant_id, zygosity}]`, `phase_sets [{variants: [variant_id…], phase}]` (trans: exactly 2), `background` (bool \| null; null = yes if configured) |
| `simulation` | optional (needed by `prepare-sim`/`simulate`): `profile` (SimulationProfile YAML), `mode` (`region` default \| `genome`), `flank_bp` (region mode, default 1500), `extra_regions_bed` (region mode, e.g. paralogous loci), `genome_contigs` (genome mode; default chr1–22, chrX, chrY), `min_alt_reads` (planting check, default 1), `mason_batch_bp` (default null: one Mason run per ploidy group; else each group is split into runs of consecutive whole segments of at most this many bp, a longer segment running alone — caps Mason's disk use, changes the reads) |
| `output` | `directory` |

**Rule keys common to all sources:** `id`, `category`, `source` (`clingen_erepo` \| `clinvar`
\| `giab` \| `explicit`), `required` (default false), `samples`, `expected {classification:
[..], reportable, level, basis?}`, `zygosity` (optional), `compound_het {phase}` (eRepo and
ClinVar only).

| `source` | categories | filters | selection mode (exactly one) |
|---|---|---|---|
| `clingen_erepo` | curated, difficult | `classification_allow` (required, strict), `genes`, `vceps`, `conditions` (MONDO ids), `held_out` (`include` default \| `exclude` \| `only`), `region_bed`, `strata`, `exclude_strata`, `variant_types` | `per_gene` (needs `genes`), `count`, `ladder_per_bin`, `erepo_ids` (Uuids); ranking: seeded hash only |
| `clinvar` | curated, difficult | `clnsig_allow` (required), `min_stars`, `clnvc`, `genes`, `conditions` (CLNDISDB ids, e.g. `MONDO:0011450`), `region_bed`, `strata`, `exclude_strata`, `variant_types` | `per_gene` (needs `genes`), `count`, `ladder_per_bin`, `variation_ids`, `rsids` (ClinVar `RS` only) |
| `giab` | standard, difficult | confident BED (containment), `region_bed`, `strata`, `exclude_strata`, `variant_types` | `counts {variant_type: n}`, `count`, `ladder_per_bin`; `level` must be `calling` |
| `explicit` | synthetic_novel | constraints only | `variants [{chrom, pos, ref, alt, gene, classification, condition, acmg_criteria, rationale}]` |

`difficult` rules need `strata` or `ladder_per_bin`. A `giab` rule whose truth VCF is the
background VCF must not target samples that carry that background (the background already is
their calling truth). Global constraints (`allowed_variant_types`,
`max_indel_length`) apply to every rule. Event length (ladder, max length): insertion
`len(alt)-1`, deletion `len(ref)-1`, delins longer allele minus a shared anchor base; ladder bins
pool insertions, deletions and delins.

## SelectionResult (`variants.selected.json`)

`schema_version`, `suite_name`, `seed`, `variants: list[TestVariant]` (sorted by reference
position), `rules: list[RuleSummary]`, `warnings`. `RuleSummary`: `rule_id`, `required`,
`samples`, `filter_counts` (ordered `[filter, candidates remaining]`), `requested` and
`selected` per bucket (gene, ladder bin, variant type, id or `total`), `rejections
[{variant_id, reason}]` (deduplication), `missing` (human-readable shortfalls), `satisfied`.
`variants.selected.tsv` has one row per TestVariant (flattened; eRepo columns `vcep`,
`erepo_classification`, `met_criteria`, `held_out`). `selection_report.md` adds a per-gene
pairing table of curated variants (positive P/LP vs negative B/LB, held-out counts, VCEP,
asymmetry flag) and notes categories below 59 variants (Roy et al. 2018).

## RuntimeConfig (YAML, `--runtime`)

How external tools are executed; independent of the suite (the same suite can run locally or in
containers). Recorded as manifest input `runtime_config`; every tool's version and container
digest are recorded in `RunManifest.tools`.

| Key | Content |
|---|---|
| `schema_version` | `1` |
| `default_runner` | `local` (default) \| `docker` |
| `docker` | `binary` (default `docker`), `extra_args` (e.g. `[--platform, linux/amd64]`) |
| `tools` | `{name: {runner, image, executable}}`; names: `bcftools`, `samtools`, `mason_simulator`, `hap.py`, `rtg`; a docker tool needs `image` (`name@sha256:…` preferred; a tag is recorded with a warning) |

## SimulationProfile (YAML, `configs/profiles/`)

| Key | Content |
|---|---|
| `schema_version`, `name` | |
| `technology` | `illumina` |
| `read_length` | bp (paired-end) |
| `fragment` | `model` (`normal` \| `uniform`), `mean`, `sd`, `min`/`max` (uniform: required) |
| `coverage` | mean depth of a diploid locus; haploid segments get half |
| `threads` | Mason threads (default 1); **part of the determinism key** (one RNG per thread) |
| `chunk_size` | Mason `--chunk-size` (≥ 65536, default 65536); also affects the output |
| `illumina` | optional Mason error-model overrides: `prob_mismatch`, `prob_mismatch_scale`, `prob_mismatch_begin`, `prob_mismatch_end`, `prob_insert`, `prob_deletion` (null = Mason default) |

## Simulation outputs (step 3)

Per sample, in the run directory:

| File | Content |
|---|---|
| `<sample>.sim/segments.tsv` | `segment, chrom, start, end, ploidy, group` — the coordinate map (1-based closed) |
| `<sample>.sim/<run>.fa(.fai)` | segment FASTA per Mason run (deleted after its run in genome mode) |
| `<sample>.sim/<run>.mason.vcf` | Mason input (constraints below); absent for a run without records |
| `<sample>.sim/background_dropped.tsv` | background records that cannot be simulated: `chrom, pos, ref, alt, gt, reason` (`overlap_on_haplotype`, `het_on_haploid_segment`, `haploid_gt_on_diploid_segment`, `crosses_segment_boundary`, `indel_at_segment_end`: mason_simulator crashes on an indel that reaches the last base of a segment, region mode) |
| `<sample>.sim/mason_<run>.log` | Mason's stderr, written while Mason runs |
| `<sample>.R1.fastq.gz`, `<sample>.R2.fastq.gz` | reads of all groups, gzip (ISA-L level 3, no file name, `mtime=0`), read names `<sample>.<run>.<n>/1\|2` |
| `<sample>.true.bam` (+ `.bai`) | Mason's true alignments lifted to genome coordinates, coordinate-sorted (`samtools sort --no-PG`); tags `oH` (haplotype, 1-based), `oP`, `oS`, `XE` (errors), `XS`, `XI`, `NM`, `MD` |
| `<sample>.simulation.json` | SimulationSpec |
| `<sample>.planting_check.json` / `.tsv` | PlantingCheck |

**Segments.** Region mode: planted variant ± `flank_bp` (clipped), plus `extra_regions_bed`,
merged per contig; genome mode: whole contigs. Male samples: chrX split at the PARs (PAR
diploid, rest haploid), chrY outside PARs haploid, chrY PARs left out (the chrX PAR copy stands for
both); female samples: no chrY. A segment shorter than the profile's longest fragment
(`max`, else `mean + 4·sd`) is an error. Background records inside segments are simulated too.

**Mason VCF** (verified against `mason_simulator` 2.0.9 / SeqAn 2.4.0 source and test runs): one sample column; FORMAT `GT` only; one ploidy per file;
`##contig` lines identical to the segment FASTA (Mason matches contigs by index); records shifted to
segment coordinates; MNV/delins **decomposed** (SNVs for the differing bases of the common prefix +
one indel anchored at the last common base; Mason silently skips records with both alleles longer
than one base); no conflicting parts on one haplotype (planted conflicts are errors, background
conflicts are dropped and logged). Checked by `validate_mason_vcf` before every run.

**Mason runs.** `<run>` is the ploidy group (`diploid`, `haploid`), or `<group>-<n>` (1-based, in
segment order) when `simulation.mason_batch_bp` splits a group into batches. Mason runs one after
another (with `TMPDIR=<sample>.sim/`); as soon as a run ends, its reads are appended to the gzip
FASTQ and its alignments lifted into the unsorted true BAM, in a background thread and in run
order while Mason simulates the next run, then its raw outputs are deleted (unless
`--keep-work`). The outputs are byte-identical to merging all groups at the end.

**Fragments and seed.** `n = ⌈coverage/2 · ploidy · L / (2 · read_length)⌉` per run (L = total
segment length incl. N); `--seed = derive_int(seed, "mason", "{sample}|{run}") mod (2^31 − 1)`.

### SimulationSpec (`<sample>.simulation.json`)

`schema_version`, `sample_id`, `mode`, `profile` (embedded SimulationProfile), `profile_sha256`,
`flank_bp` (null in genome mode), `segments` (path), `runs: list[MasonRun]`, `decomposed`
(planted variant ids split for Mason), `background_dropped`, `background_dropped_log`,
`fastq_r1`, `fastq_r2`, `true_bam`, `planting_check`, `planting_complete`.
`MasonRun`: `group` (`diploid`\|`haploid`), `name` (batch name `<group>-<n>`; null without
batches), `segment_names` (segments of the batch; null without batches), `ploidy`, `segments`,
`length`, `fragments`, `seed`, `vcf`, `vcf_records`, `args` (complete `mason_simulator` arguments,
paths relative to its cwd).

### PlantingCheck (`<sample>.planting_check.json`)

`schema_version`, `sample_id`, `min_alt_reads`, `window_pad`, `complete` (all `ok`),
`variants: list[VariantPlanting]` sorted by id. `VariantPlanting`: `variant_id`, `chrom`, `pos`,
`ref`, `alt`, `zygosity`, `carrying_haplotypes` (Mason haplotypes: 1 ⇔ GT `1|0`, 2 ⇔ `0|1`,
hemizygous: 1), `expected_vaf` (het 0.5, hom/hemi 1.0), `spanning_reads`, `alt_reads`,
`ref_reads`, `other_reads`, `vaf`, `by_haplotype [{haplotype, reads, alt, ref, other}]`,
`status`.

Window: REF span extended to the right by the indel's repeat shift, padded on both sides by
`window_pad + max(len(REF), len(ALT))`. Only primary reads spanning the whole window count. A read
whose sequence contains the whole ALT window sequence but not the REF one (or the reverse) is
classified by that (for a tandem duplication inside a repeat Mason's CIGAR can scatter the
inserted bases over many 1–2 bp gaps, so the aligned bases no longer match either allele);
otherwise the read bases aligned inside the window are compared with the REF and ALT window
sequences (exact, then edit distance; ties = other). `clinval recheck-planting` repeats the check
from an existing true BAM (no Mason run) and rewrites the check files and `planting_complete`. Indels longer than 20 bp (beyond the reference implementation) have a window longer than a read, and indels that can shift in a repeat have a wider window that few of Mason's reads span (their starts are shifted); for both, a read that overlaps the window without spanning it is classified by marker k-mers (40 bp) of the local haplotypes (20 bp of flank, the allele, the reference over which the indel can shift in a repeat, and 20 bp after it), keeping only k-mers absent from the other allele's local haplotype: ALT if the read contains an ALT marker and no REF marker, and vice versa. The length-dependent padding is needed because Mason realigns reads with
edit-distance scoring and free end gaps: an indel closer to a read end than its own length is
represented by a shifted read start instead of a gap (observed for a 3 bp insertion), so such
reads are uninformative. `status`: `ok` if every carrying haplotype has ≥ `min_alt_reads` ALT
reads, `no_reads` if a carrying haplotype has no spanning read, else `missing_on_haplotype`.

## CanonicalFinding (tertiary output, normalised)

| Field | Type |
|---|---|
| `sample_id`, `chrom`, `pos`, `ref`, `alt` | |
| `gene`, `hgvs_c`, `hgvs_p` | str \| null |
| `classification` | Classification \| null |
| `reported` | bool |
| `acmg_criteria` | list[str] |
| `tool`, `tool_version` | str |

## OutcomeRecord and OutcomeSet (`sut/<pipeline>/compare/<sample>.outcomes.json`)

`OutcomeRecord`: `variant_id` (null for a false report), `sample_id`, `pipeline`, `chrom`, `pos`,
`ref`, `alt`, `evaluation_level` (null for a false report), `outcome_class`, `stage` (last check
passed: `simulated`, `covered`, `mapped`, `called`, `passed`, `genotyped`, `annotated`,
`classified`, `reported`; null if none), `expected_classification`, `expected_reportable`,
`observed_classification`, `observed_reported`, `gene` (of the matched finding), `evidence`
(dict of strings, numbers, booleans, lists).

Evidence of a planted variant: `planting_status`, `true_alt_reads`, `aligned_alt_reads`,
`aligned_spanning`, `aligned_alt_mapq0` (null without an alignment), `happy_decision` (`BD` of the
truth allele: TP, FN; null if hap.py has no record), `happy_kind` (`BK`: gm, am, lm, .),
`happy_match` (`exact` or `window` fallback), `nonpass_filters` (FILTER of non-PASS raw-VCF
records of the allele), `calls_in_window`, `finding_found`, `acmg_criteria`. Evidence of a false
report: `call` (`fp_call`, `true_variant` from the exact hap.py match; otherwise
`outside_scored_regions` = not in the scored segments and planted windows, or `no_exact_match` =
scored but no hap.py record with exactly this allele), `happy_decision`, `planted_negative`,
`planted_variant_id`, `acmg_criteria`. A reported finding inside the window of a planted variant
that hap.py found but whose allele was not annotated is the same failure in another
representation: it is not a false report and is listed in that variant's evidence as
`split_report` (`chrom:pos:ref:alt`, comma-separated). A reported planted variant evaluated at `calling` level
only is not a false report (it has no expected report decision).

`OutcomeSet`: `schema_version`, `sample_id`, `pipeline`, `parameters` (`window_pad`,
`mismapped_fraction`, `edge_bp`), `records` (one per planted variant, by `variant_id`),
`false_reports` (sorted by position). `outcomes.tsv` lists the records of all compared samples.

`OutcomeClass`, ordered along the pipeline (first failing stage wins):

| Class | Meaning |
|---|---|
| `not_simulated` | variant not present in simulated reads (true alignments) — framework error |
| `not_covered` | insufficient simulated depth at locus |
| `mismapped` | allele not called and reads carrying it not aligned to the true locus (e.g. paralogues); a called allele is never `mismapped`, because aligners clip reads around long indels |
| `not_called` | absent from raw VCF |
| `filtered` | in raw VCF but removed / FILTER≠PASS |
| `genotype_mismatch` | called with wrong zygosity |
| `not_annotated` | lost or wrongly annotated during annotation |
| `misclassified` | classification outside `acceptable_classifications` |
| `report_mismatch` | correct classification but `reported` ≠ `expected_reportable` (also when the system under test reports by a different policy, e.g. recessive carriers) |
| `correct` | all scored layers correct |
| `false_report` | (negative controls / off-target) reported although it must not be |

## AcceptanceConfig, MetricsReport, DiffReport

`AcceptanceConfig` (YAML, `configs/acceptance/`, separate from the suite so that it does not
invalidate a run directory): `schema_version`, `name`, `criteria: [Criterion]`. `Criterion`:
`metric`, `subset` (dict, e.g. `{held_out: "true"}`), and `min_point`, `min_lower` (lower Wilson
bound) for proportions or `max_count` for counts. A proportion criterion uses the row whose subset
equals its own; a count criterion sums every row whose subset contains its own.

`MetricsReport` (`metrics.json`): `schema_version`, `pipeline`, `runs`, `proportions: [MetricRow
{metric, subset, value: Proportion}]`, `counts: [CountRow {metric, subset, count}]`,
`acceptance`, `criteria: [CriterionResult {criterion, status: pass|fail|not_evaluable,
observed}]`. Metrics and their subsets:

| Metric | Over | Subsets |
|---|---|---|
| `report_sensitivity` | planted, `expected_reportable` | all, `held_out`, `suite` |
| `report_specificity` | planted, not reportable (B/LB, VUS, recessive carriers) | all, `group`, `held_out`, `group`+`held_out`, `suite` |
| `class_exact`, `class_acceptable`, `class_tier` | planted with an observed class | all, `group` |
| `criterion_recall`, `criterion_precision` | ACMG codes without strength, planted with a finding | all, `criterion` |
| `planted_calling` | planted, not failed before `genotyped` | all, `variant_type`, `suite` |
| `calling_recall`, `calling_precision` | hap.py k/n summed over samples | `type`+`stratum`; recall also `contig_group` |
| `stage_counts` (count) | planted | `outcome_class`+`group` |
| `false_reports` (count) | false reports | `call`+`planted_negative`, `suite`+`sample_id` |
| `false_positive_calls` (count) | hap.py query FP | `contig_group` (autosomes, chrX_PAR, chrX_nonPAR, …; PAR split with `--par-bed`) |

`group`: P/LP (reportable), recessive_carrier (P/LP, not reportable), B/LB, VUS, calling_only.
`held_out`: from `held_out.tsv`; held out but approved before the reference snapshot =
`unknown`, no entry = `unknown`. A variant lost before the report counts as not reported.

`DiffReport` (`diff.json`): `schema_version`, `baseline`, `candidate`, `noise`, `runs`,
`changes: [Change {kind, run, sample_id, variant_id, chrom, pos, ref, alt, gene, baseline,
candidate, noise}]` (kind: `new_failure`, `recovered`, `stage_change`, `class_change`,
`false_report_added`, `false_report_removed`; baseline/candidate = `outcome/class/reported`),
`calling: [CallingDelta {run, sample_id, contig_group, counts}]` (truth TP/FN and query FP of
both), `concordance: [Concordance {run, sample_id, contig_group, identical, only_baseline,
only_candidate}]` (PASS calls; empty when both pipelines share the secondary outputs), `verdict`
(`regression` if a `new_failure` or `false_report_added` is not noise), `acceptance`, `criteria`.

## FidelityConfig (YAML, `configs/fidelity/`)

Separate from the suite; paths relative to the file; unknown keys are errors.
`schema_version`, `name`, `reference` (FASTA used by the pipeline), `confident_bed` (regions of
hap.py, e.g. GIAB v4.2.1 benchmark BED), `qc_regions` (regions of `samtools stats` /
`bcftools stats`; null = `confident_bed`), `stratifications` (stratum name → BED, plain or .gz),
`datasets: list[FidelityDataset]` (≥ 2, unique `id`), `baseline` (dataset id, the real data),
`recall_gap` (flag threshold, null = no flags), `happy {engine: vcfeval|xcmp, threads, pass_only,
vcfeval_template, gender}` (`gender`: `male`|`female`|null, hap.py `--gender`; null = hap.py
guesses and treats chrX as diploid), `planting` (run directories), `threads`.
`FidelityDataset`: `id`, `label`, `alignment` (BAM/CRAM), `vcf`, `truth_vcf`.

## FidelityResult (`fidelity.json`)

`schema_version`, `name`, `baseline`, `datasets {id: label}`, `pass_only`, `recall_gap`,
`planting` (PlantingSummary \| null), `read_qc: list[ReadQc]`, `callset_qc: list[CallsetQc]`,
`calling: list[CallingMetric]`, `gaps: list[RecallGap]`, `missing` (raw outputs not found).

| Model | Fields |
|---|---|
| `Proportion` | `k`, `n`, `value`, `lower`, `upper` (95 % Wilson; null when `n` = 0) |
| `ReadQc` | `dataset`, `reads_total`, `reads_mapped`, `reads_properly_paired`, `reads_mq0` (fractions of reads), `insert_size_mean`, `insert_size_sd`, `error_rate` (mismatches / bases mapped), `indel_rate` (indel events from `IC` / bases mapped), `coverage_mean`, `coverage_ge20`, `fold80` (mean / 20th percentile of depth) |
| `CallsetQc` | `dataset`, `snvs`, `indels`, `ti_tv`, `het_hom` |
| `CallingMetric` | `dataset`, `type` (SNP\|INDEL), `stratum` (`*` = all), `subtype` (hap.py, `*` = all), `recall` (TRUTH.TP of TP+FN), `precision` (QUERY.TP of TP+FP), `f1` |
| `RecallGap` | `dataset`, `type`, `stratum`, `recall`, `baseline_recall`, `gap` (recall − baseline), `flagged` (\|gap\| > `recall_gap`) |
| `PlantingSummary` | `runs`, `samples`, `variants`, `complete` (Proportion of `ok`), `by_status`, `by_zygosity [{zygosity, variants, vaf_median, vaf_min, vaf_max, alt_reads_median, alt_reads_min}]` |

Error and indel rates are measured against the reference and include the genome's true variants:
they compare datasets of the same genome, they are not sequencing error rates. Coverage is
computed from the `COV` histogram over the length of `qc_regions` (uncovered bases count as 0;
depths above 1000 count as 1000). `planting.tsv`: one row per planted variant and sample
(`run`, `sample_id`, `variant_id`, `zygosity`, `expected_vaf`, `spanning_reads`, `alt_reads`,
`vaf`, `status`).

## RunManifest

`schema_version`, `run_id`, `created_at`, `framework {version, git_commit, dirty, python_packages}`,
`config {path, sha256}`, `seed`, `inputs [{role, path, sha256}]`, `outputs [{role, path, sha256}]`,
`tools [{name, version, container_digest}]`, `databases [{name, release, sha256}]`,
`steps [{name, started_at, wall_s, peak_rss_mb, exit_code}]`, `warnings [str]`.

Written as `run_manifest.json` in the run directory and extended by every command writing into
it (a manifest from another config hash or seed is refused). Output paths are relative to the run
directory. `run_id` = UTC timestamp + first 8 hex digits of the config hash. Framework steps
report the process-lifetime peak RSS; tool steps (e.g. `bcftools norm (positive)`) the child's;
for tools run in docker `peak_rss_mb` is null (the measured process is the docker client) and
`tools[].container_digest` holds the image digest (pinned) or local image id.
`git_commit`/`dirty` are null outside a git checkout or before the first commit.
`python_packages` maps the Python libraries whose version changes output bytes to their installed
version: `pysam` (htslib: BAM, VCF, BGZF) and `isal` (ISA-L gzip of the FASTQ). A command that
extends an existing manifest adds missing entries and warns when a recorded version differs.
