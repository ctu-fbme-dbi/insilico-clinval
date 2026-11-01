# User manual

How to install the framework, get the data, run every step and read the outputs. For a
one-minute end-to-end run on bundled data, see the quickstart in the [README](../README.md).

```
suite.yaml ──► 1 select ──► 2 build-truth ──► 3 simulate ──► 4-5 run-sut ──► 6 compare ──► metrics / diff / report
               variants      truth VCF +        Mason reads,    pipeline under    outcomes,
                             expected outcomes  true BAM, check test (adapters)   failure stages
```

---

## 1. Prerequisites

| Requirement | Needed for | Check |
|---|---|---|
| Python ≥ 3.10 | everything | `python3 --version` |
| `bcftools` | step 2 (`build-truth`) | `bcftools --version` |
| `samtools` | step 3 (`simulate`) | `samtools --version` |
| Docker | Mason (step 3), hap.py (step 6); optional for bcftools/samtools | `docker info` |
| ~10 GB disk | reference data for the demo suite (GRCh38 is 3.1 GB uncompressed) | |
| Nextflow | only for nf-core/sarek (step 4) | `nextflow -version` |

On macOS: `brew install bcftools samtools`, and start Docker Desktop (`open -a Docker`) before
running step 3. Instead of Docker, the tools can run from PATH (see section 8).

## 2. Installation

```bash
cd insilico-clinval/          # the repository root
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
clinval --help
```

Activate the environment (`source .venv/bin/activate`) in every new terminal.

## 3. Build the tool images (once)

```bash
docker build -t clinval-mason:2.4.0 containers/mason                              # step 3 (Mason)
docker build --platform linux/amd64 -t clinval-happy:0.3.15 containers/happy       # step 6 (hap.py)
docker build -t clinval-nirvana:3.18.1 containers/nirvana                          # BIAS-2015 with Nirvana
docker build -t clinval-bias2015:v3.0.0 containers/bias2015                        # BIAS-2015
docker build -f containers/clinval/Dockerfile -t clinval:dev .                     # optional: bcftools/samtools in docker
```

`clinval-mason` compiles `mason_simulator` from the SeqAn 2.4.0 sources (about a minute). The
quickstart needs only the first two; see `containers/README.md` for all images.

## 4. Quick start on the bundled fixtures (no downloads)

The repository contains a tiny synthetic genome and test data (`tests/fixtures/`).
`examples/quickstart/run_quickstart.sh` runs all six steps on it with a toy pipeline (see the
README). Steps 1–3 alone take seconds and are the best way to check the installation:

```bash
# steps 1+2: select variants and write the truth
clinval build-truth --suite tests/fixtures/suite_erepo.yaml --out runs/fixture_erepo

# step 3: simulate reads (Mason in docker)
clinval simulate --suite tests/fixtures/suite_erepo.yaml --out runs/fixture_erepo \
                 --runtime configs/runtime/docker-local-build.yaml

# look at the result
cat runs/fixture_erepo/selection_report.md
cat runs/fixture_erepo/positive.planting_check.tsv
```

Every line of `positive.planting_check.tsv` should end with `ok`.

## 5. Download the reference data

All downloads go to `data/` (git-ignored) and are verified with SHA-256.

```bash
python scripts/fetch_data.py --list          # what is available, pinned or not

python scripts/fetch_data.py \
    --only grch38,giab_hg002_v4.2.1,clingen_erepo \
    --allow-unpinned
```

| Dataset | Content | Size |
|---|---|---|
| `grch38` | GRCh38 no-alt analysis set → `data/reference/GRCh38_no_alt_analysis_set.fa` (+ `.fai`) | 0.9 GB download, 3.1 GB unpacked |
| `giab_hg002_v4.2.1` | HG002 benchmark VCF + BED (background genome) | 170 MB |
| `clingen_erepo` | ClinGen eRepo export → `data/clingen/erepo_<YYYY-MM-DD>.tsv` | 30 MB |
| `clinvar_20260928` | ClinVar VCF (only for ClinVar-based suites) | 200 MB |
| `giab_cmrg_v1.00`, `giab_strata_v3.3` | difficult-region BEDs (only for `difficult` rules) | small |
| `bias2015_erepo_benchmark_v2.2.0` | held-out reference for BIAS-2015: the eRepo V2.2.0 benchmark variants (Eisenhart et al. 2025, Additional file 2), lifted to GRCh38 with CrossMap (docker) | 4 MB |
| `bias2015_hg38` (+ `bias2015_hg38_vep`) | BIAS-2015 v3.0.0 data bundle 2026.03.01 (+ VEP-specific files) | 0.8 GB |
| `nirvana_grch38`, `nirvana_grch38_updates` | Nirvana 3.18.1 GRCh38 data (Illumina Downloader; build `clinval-nirvana` first) + Bitscopic's ClinVar 2026-02 and AlphaMissense; the stale bundled ClinVar files are removed | tens of GB |
| `vep_cache_grch38_115`, `vep_hg38_annotations`, `gnomad_v4.1_exomes` | for BIAS-2015 in VEP mode: VEP 115 cache + FASTA, ClinVar 2026-02-01, REVEL, AlphaMissense, gnomAD v4.1 exomes reduced to AF/AF_joint/AC/AN | ~26 GB + 1.3 GB + ~300 GB download |

Notes:

- `--allow-unpinned` is needed for the eRepo export: it changes daily, so it has no fixed
  checksum. Its hash is still recorded in `data/fetch_report.json` and in every run manifest.
- Re-running the command does not download again; existing files are only verified.
- A checksum mismatch stops with an error; the bad file is renamed `*.sha256-mismatch`.
- **Error "self signed certificate in certificate chain"** (TLS-inspecting proxy): export the
  system certificates and pass them:

  ```bash
  security find-certificate -a -p /Library/Keychains/System.keychain \
      /System/Library/Keychains/SystemRootCertificates.keychain > /tmp/ca.pem   # macOS
  python scripts/fetch_data.py --ca-bundle /tmp/ca.pem --only ... --allow-unpinned
  ```

## 6. Run the demo suite on real data

`configs/suites/demo.yaml`: ClinGen eRepo pathogenic/likely pathogenic variants (sample
`positive`) and benign/likely benign variants (sample `negative`) on CDH1, TP53, MYH7, PTEN, LDLR
and MECP2, planted into the HG002 genome, simulated at 30x around each variant.

**6.1 Point the suite to your eRepo download.** Edit `configs/suites/demo.yaml`:

```yaml
sources:
  clingen_erepo:
    file: ../../data/clingen/erepo_2026-10-07.tsv   # <- your file name
    release: "2026-10-07"                           # <- the same date
```

**6.2 Validate the suite** (schema, references between rules and samples, input files exist):

```bash
clinval validate-suite configs/suites/demo.yaml
```

**6.3 Step 1 — select variants** (optional on its own; `build-truth` includes it):

```bash
clinval select --suite configs/suites/demo.yaml
less runs/demo/selection_report.md
```

**6.4 Steps 1+2 — truth set** (about 2 minutes, dominated by merging the HG002 background):

```bash
clinval build-truth --suite configs/suites/demo.yaml
```

**6.5 Step 3 — simulation** (a few seconds in region mode):

```bash
clinval simulate --suite configs/suites/demo.yaml \
                 --runtime configs/runtime/docker-local-build.yaml
```

Without `--out`, results go to `output.directory` of the suite (`runs/demo`). Use
`--sample positive` (repeatable) to simulate only some samples, and `--keep-work` to keep Mason's
raw outputs. `clinval prepare-sim` runs only step 3a (writes the Mason inputs without simulating),
which is useful for inspecting them.

## 7. What you get (run directory)

### Step 1 — selection

| File | Content |
|---|---|
| `selection_report.md` | per rule: candidates left after each filter, selected per gene, rejections, shortfalls; per-gene table of positive vs negative variants; skipped eRepo records by reason |
| `variants.selected.json` / `.tsv` | every selected variant with its source (database, release, record id), classification, VCEP, ACMG criteria, `held_out` |

### Step 2 — truth (per sample)

| File | Content |
|---|---|
| `<sample>.truth.vcf.gz` (+ `.tbi`) | planted variants + HG002 background, phased genotypes; this is the calling truth |
| `<sample>.expected_outcomes.json` | per planted variant: expected classification, acceptable classifications, whether it must be reported |
| `<sample>.sample.json` | sex, zygosity/haplotype of each variant, background summary |
| `<sample>.background_removed.tsv` | HG002 records removed because they were too close to a planted variant |

### Step 3 — simulation (per sample)

| File | Content |
|---|---|
| `<sample>.R1.fastq.gz`, `<sample>.R2.fastq.gz` | **the input for the pipeline under test** |
| `<sample>.true.bam` (+ `.bai`) | where every read really comes from (Mason's true alignments, genome coordinates; tag `oH` = haplotype) |
| `<sample>.planting_check.tsv` / `.json` | per planted variant: reads spanning it, reads carrying it per haplotype, VAF, status |
| `<sample>.simulation.json` | the simulation profile, segments and every Mason run (arguments, seed, number of fragments) |
| `<sample>.sim/` | `segments.tsv` (simulated regions), segment FASTA, Mason input VCFs, Mason logs, `background_dropped.tsv` |

### Steps 4–6 — system under test and comparison

After `run-sut` (a pipeline config from `configs/adapters/`), compare its outputs with the truth:

```bash
clinval compare --suite configs/suites/<suite>.yaml --out runs/<suite> \
    --pipeline configs/adapters/<pipeline>.yaml --runtime configs/runtime/docker.yaml \
    --pipeline-reference <FASTA the pipeline aligned to>   # decodes CRAM; default: the suite's
```

Outputs in `runs/<suite>/sut/<pipeline>/compare/`:

| File | Content |
|---|---|
| `<sample>.outcomes.json` | one record per planted variant (outcome class, last stage passed, evidence) and every false report |
| `outcomes.tsv` | the same for all compared samples, one row each |
| `<sample>.calling.json` | hap.py recall and precision (PASS calls) per variant type and stratum: `planted`, `giab_confident` (the suite's background BED or `--confident-bed`) and the suite's `sources.strata` |
| `<sample>.regions.bed` | where calls were scored: the simulated segments (ends trimmed by `--edge-bp`) plus a window around every planted variant; the truth is complete there by construction |
| `<sample>.happy/` | raw hap.py output |

Outcome classes and their order are in `docs/schemas.md`. A false report's `call` says whether
the reported allele was a false-positive call (`fp_call`) or a real variant (`true_variant`); a
report without an exact allele match in the hap.py output is `outside_scored_regions` (not in the
simulated segments or planted windows) or `no_exact_match` (scored, but hap.py has no record with
exactly this allele, e.g. a long insertion represented differently). A reported part of a planted
variant that was called in another representation (a delins called as an insertion and an SNV)
is not a false report: the planted variant's record carries it as `split_report`. On chrX, where
HG002 v4.2.1 has no background, every call except a planted variant is a false positive.

### Metrics, acceptance criteria and regression diff

```bash
clinval metrics --pipeline configs/adapters/<pipeline>.yaml --run runs/<suite1> --run runs/<suite2> \
    --out results/<pipeline> --acceptance configs/acceptance/article-default.yaml \
    --par-bed configs/regions/grch38_par.bed
clinval diff --baseline configs/adapters/<baseline>.yaml --candidate configs/adapters/<candidate>.yaml \
    --run runs/<suite> --out results/diff --par-bed configs/regions/grch38_par.bed
clinval report --metrics results/<pipeline> --diff results/diff --out results/report.md
```

`metrics` writes `metrics.json` and `metrics.tsv` (every proportion with k, n and the 95 % Wilson
interval; counts); `diff` writes `diff.json`, `diff.tsv` (one row per change) and `diff.md`.
Metric definitions are in `docs/schemas.md`. `--noise` names a second run of the baseline;
changes also seen there are marked as noise and do not make a regression.

### Held-out status of an existing run

```bash
clinval held-out --suite configs/suites/<suite>.yaml --out runs/<suite> \
    --reference data/clingen/bias2015_erepo_v2.2.0_benchmark.grch38.vcf --release v2.2.0 \
    --approved-after 2025-02-02
```

Writes `held_out.tsv` (one row per eRepo variant). Run it before `clinval compare`, which copies
the flag into its outcomes.

### Always

| File | Content |
|---|---|
| `run_manifest.json` | provenance: framework version and git commit, config hash, seed, SHA-256 of all inputs and outputs, tool versions and container digests, time of every step, warnings |
| `clinval.log` | detailed log |

### Reading the planting check

```
variant_id                     zygosity expected_vaf vaf    spanning_reads alt_reads ... alt_by_haplotype status
curated:chr17:7674200:TGATGGTGAG:T het  0.5          0.5    28             14        ... 1:0/14,2:14/14   ok
curated:chrX:154031329:G:A     hemi     1.0          1.0    17             17        ... 1:17/17          ok
```

`alt_by_haplotype` `2:14/14` means: 14 reads from haplotype 2 span the variant and all 14 carry
it. A heterozygous variant should be on one haplotype only, a homozygous one on both, a
hemizygous one (chrX in a male) on haplotype 1 at about half depth. Status: `ok`;
`missing_on_haplotype` (reads exist but do not carry the variant); `no_reads` (no read spans it,
e.g. too low coverage). Any status other than `ok` makes `simulate` exit with code 4.

Useful commands:

```bash
samtools view runs/demo/positive.true.bam chr17:7674150-7674250 | head
zcat runs/demo/positive.R1.fastq.gz | head -8          # gzcat on macOS
python -m json.tool runs/demo/run_manifest.json | less
```

## 8. Configuration

### 8.1 The suite (`configs/suites/*.yaml`)

The suite defines everything that affects the results; the same suite + seed gives the same
outputs. Main sections (full reference: `docs/schemas.md`, *SuiteConfig*; annotated example:
`configs/suites/example.yaml`):

```yaml
schema_version: 1
suite: { name: my-suite, version: "0.1", seed: 12345 }     # seed drives every random choice
reference: { build: GRCh38, fasta: ../../data/reference/GRCh38_no_alt_analysis_set.fa }
sources:                     # where variants come from
  clingen_erepo: { file: ..., release: "..." }
constraints: { min_distance_bp: 500, allowed_variant_types: [snv, insertion, deletion, delins], max_indel_length: 20 }
selection:                   # rules, processed in order
  - id: erepo_plp
    category: curated
    source: clingen_erepo
    classification_allow: [pathogenic, likely_pathogenic]
    genes: [TP53, PTEN]
    per_gene: 2
    required: true
    expected: { classification: [pathogenic, likely_pathogenic], reportable: true, level: both }
    samples: [positive]
background: { name: GIAB_HG002, vcf: ..., bed: ..., release: v4.2.1, sex: male }
samples:
  - { sample_id: positive, role: positive, sex: male }
simulation:
  profile: ../profiles/illumina_2x150_30x.yaml
  mode: region              # region: around planted variants; genome: whole chromosomes
  flank_bp: 1500            # region mode: bases simulated on each side of a variant
  # extra_regions_bed: ...  # region mode: also simulate e.g. paralogous loci
output: { directory: ../../runs/my-suite }
```

Relative paths are resolved from the suite file's directory. Unknown keys are errors.

Variant sources and selection modes:

| `source` | typical use | select by |
|---|---|---|
| `clingen_erepo` | expert-curated P/LP and B/LB (reference choice) | `per_gene`, `count`, `ladder_per_bin`, `erepo_ids` |
| `clinvar` | ClinVar, filtered by review stars | `per_gene`, `count`, `ladder_per_bin`, `variation_ids`, `rsids` |
| `giab` | calling-only variants from a GIAB truth set | `counts`, `count`, `ladder_per_bin` |
| `explicit` | your own synthetic novel variants | `variants: [...]` |

Per sample you can set `default_zygosity` (`het`/`hom`), per-variant `overrides`, and
compound-heterozygous `phase_sets` (`cis`/`trans`).

### 8.2 Simulation profiles (`configs/profiles/`)

| Profile | Reads | Fragments | Depth |
|---|---|---|---|
| `illumina_2x150_30x.yaml` | 2 × 150 bp | 400 ± 60 bp | 30x |
| `illumina_2x150_30x_t8.yaml` | 2 × 150 bp | 400 ± 60 bp | 30x, 8 Mason threads (the article's suites) |
| `illumina_2x150_100x_panel.yaml` | 2 × 150 bp | 400 ± 60 bp | 100x |

Optional keys: `threads` (Mason threads; changing it changes the reads), `illumina` error-model
parameters. Reference: `docs/schemas.md`, *SimulationProfile*.

**Region vs genome mode.** Region mode simulates only the planted loci (± flank): fast, small
FASTQ, good for development and change control, but the aligner never sees reads from elsewhere
in the genome (add paralogues via `extra_regions_bed`). Genome mode simulates whole chromosomes:
realistic, but a 30x genome is ~100 GB of FASTQ and hours of compute (see
`docs/acmg_sf32_server.md` for measured disk and time).

### 8.3 Runtime (`--runtime`, `configs/runtime/`)

Where external tools run. Does not change the suite; recorded in the manifest.

| File | bcftools / samtools | mason_simulator, hap.py, Nirvana, VEP, BIAS-2015 |
|---|---|---|
| (none) or `local.yaml` | PATH | PATH (must be installed) |
| `docker-local-build.yaml` | PATH | docker images built from `containers/` |
| `docker.yaml` | docker `clinval:dev` | docker images built from `containers/` |

Nextflow (Sarek) always runs on the host. Every command that runs external tools accepts
`--runtime`. Tag-only images (as here) produce a manifest
warning; for article runs, push the images to a registry and pin them as `name@sha256:…`.

### 8.4 Pipeline configs (`--pipeline`, `configs/adapters/`)

The system under test, kept separate from the suite so that several pipelines can be run on the
same simulated data in one run directory (outputs go to `runs/<suite>/sut/<pipeline name>/`).
Schema: `docs/json_schema/PipelineConfig.schema.json`; check a config with
`clinval validate-pipeline FILE`.

```yaml
schema_version: 1
name: my-pipeline                 # names the output directory; must be unique per run directory
secondary:                        # FASTQ -> alignment + VCF
  adapter: command                # sarek | command | precomputed
  command:
    command:
      tool_name: my-caller
      version: "1.2"              # recorded in the manifest (black boxes have no --version)
      runner: docker              # or local
      image: my-caller@sha256:...
      argv: [run_caller, "{fastq_r1}", "{fastq_r2}", "{fasta}", "{outdir}", "{sample}"]
    outputs: { raw_vcf: calls.vcf.gz, alignment: aligned.cram }   # relative to {outdir}
tertiary:                         # VCF -> classified, reported findings
  adapter: command                # bias2015 | command | none
  input: filtered                 # or raw
  pass_only: true
  report_classifications: [pathogenic, likely_pathogenic]   # what counts as "reported"
  command:
    command: { tool_name: my-acmg, version: "2.0", runner: local,
               argv: [my_acmg, "{input_vcf}", "{output}"] }
    output: "{sample}.acmg.tsv"
    parser:                       # maps the tool's table to the framework's findings
      format: tsv                 # tsv | csv | json
      columns: { chrom: CHROM, pos: POS, ref: REF, alt: ALT, classification: CLASS }
      classification_map: { Pathogenic: pathogenic, VUS: uncertain_significance, Benign: benign }
```

Placeholders: secondary `{sample} {fastq_r1} {fastq_r2} {fasta} {outdir}`; tertiary
`{input_vcf} {output} {sample} {fasta} {workdir}`. Paths in the config are relative to the config
file. Working examples: `examples/quickstart/toy_pipeline.yaml.in` (both stages as commands),
`configs/adapters/example_blackbox_classifier.yaml` (Sarek + a closed classifier) and the
article's Sarek + BIAS-2015 configs.

`run-sut --stage secondary|tertiary` runs one stage; `--sarek-results DIR` adopts a Sarek run made
outside the framework (release, parameters and FASTQ files are checked against the config).

## 9. Exit codes

| Code | Meaning | Typical cause |
|---|---|---|
| 0 | success | |
| 1 | invalid input or tool failure | missing file, invalid YAML, Docker not running, tool error |
| 2 | selection unsatisfied | a `required` rule found too few variants, or a sample is empty (see `selection_report.md`) |
| 3 | truth encoding failed | REF base does not match the reference FASTA |
| 4 | simulation failed | planted variants overlap on a haplotype, or the planting check is incomplete |
| 5 | system under test failed | Sarek, the annotator or the classifier exited with an error (`run-sut`) |
| 6 | comparison failed | hap.py exited with an error or wrote unreadable output (`compare`) |
| 7 | rejected | an acceptance criterion failed (`metrics --fail-on-reject`) or a regression was found (`diff --fail-on-regression`) |

Partial results are always written before a non-zero exit.

## 10. Troubleshooting

| Message | Fix |
|---|---|
| `failed to connect to the docker API` | start Docker Desktop (`open -a Docker`) and retry |
| `... was written for a different suite config or seed; use a fresh output directory` | the run directory belongs to another version of the suite; use a new `--out` or delete the old directory |
| `... not found (run clinval build-truth first)` | run `build-truth` with the same suite and `--out` before `simulate` |
| `the suite has no simulation section` | add a `simulation:` block (section 8.1) |
| `missing required inputs` | download the data (section 5) or fix paths in the suite |
| `segments shorter than the longest fragment` | increase `simulation.flank_bp` |
| `self signed certificate in certificate chain` | `--ca-bundle` (section 5) |
| selection shortfall (exit 2) | fewer variants available than requested: lower `per_gene`, add genes, or relax filters |

Step timings in `run_manifest.json` do not include time when the computer was asleep.

## 11. Run the tests

```bash
pytest                                                    # tests needing Mason or hap.py are skipped
CLINVAL_TEST_MASON_IMAGE=clinval-mason:2.4.0 \
CLINVAL_TEST_HAPPY_IMAGE=clinval-happy:0.3.15 pytest      # all tests, Mason and hap.py in docker
CLINVAL_TEST_SAREK=1 pytest tests/integration/test_sarek.py   # opt-in: nf-core/sarek (Nextflow)
ruff check . && ruff format --check . && mypy             # code quality
```
