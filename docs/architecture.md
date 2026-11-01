# Architecture

Thin Python orchestration and evaluation layer over established tools. The framework does **not**
implement its own read simulator, aligner, variant caller or ACMG classifier; it wraps them.
Own logic lives where the article's contributions are: variant selection, the expected-outcome
model, finding-level comparison, failure-stage attribution, regression diff and provenance.

## Data flow

```
suite.yaml
   │
   ▼
[variants] ── candidate TestVariants (with provenance) ─┐
   │                                                     │
   ▼                                                     │
[truth] ── truth VCF (phased GT) + expected_outcomes.json + sample spec
   │
   ▼
[simprep] ── Mason-compatible VCF (+ optional mini reference & coordinate map for region mode)
   │
   ▼
[simulate] ── mason_simulator ─► R1/R2 FASTQ + true read alignments (SAM/BAM) + checksums
   │
   ▼
[adapters: secondary] ── BAM, raw VCF, filtered VCF      (system under test, black box)
   │
   ▼
[adapters: tertiary] ── annotated VCF, interpretation report → CanonicalFinding[]
   │
   ▼
[compare] ── VCF-level (hap.py / vcfeval, stratified) + finding-level OutcomeRecord[]
   │                                  with failure-stage attribution
   ▼
[metrics] ── sensitivity (P/LP), specificity (B/LB), CIs, acceptance criteria
[diff]    ── run A vs run B (regression)
[report]  ── JSON + HTML/Markdown + paper figures
```

Every step appends to the **run manifest** (`schemas.md` → RunManifest).

## Modules

| Module | Step | Responsibility |
|---|---|---|
| `variants` | 1 | Candidate selection per category: `standard` (GIAB truth subset), `curated` (ClinGen eRepo — reference choice — or ClinVar), `difficult` (stratification BEDs such as CMRG, indel length ladder), `synthetic_novel` (explicit list; generator later). eRepo: GRCh38 mapping from genomic HGVS, allele merging, `held_out` against a reference release. Deterministic ranking, overlap/min-distance deduplication. |
| `truth` | 2 | Zygosity and phase assignment (incl. compound heterozygotes in cis/trans, sex-aware chrX/chrY), REF check against FASTA, background merge (GIAB HG002 v4.2.1, or the assembly-based v5.0q with the opt-in record filters), normalisation, truth VCF writing, expected-outcome metadata. |
| `simprep` | 3a | Segments (region mode: variants ± flank + extra BED, e.g. paralogues; genome mode: whole contigs; male chrX/chrY split at the PARs) per ploidy group; segment FASTA; Mason VCF (GT only, shifted, MNV/delins decomposed, per-haplotype conflict checks, background records that cannot be represented dropped and logged); `validate_mason_vcf`; coordinate map `segments.tsv`. |
| `simulate` | 3b | One `mason_simulator` run per sample and ploidy group (profile → arguments, fragment count, derived seed); FASTQ merge (deterministic gzip); true alignments lifted to genome coordinates, sorted and indexed; planting-completeness check per haplotype from the true alignments (read sequence first, then the aligned window bases; `recheck-planting` repeats it from an existing true BAM). |
| `adapters` | 4–5 | Black-box contract. Generic YAML-configured command adapter (container, command, declared outputs) + reference adapters for the open-source tools evaluated in the paper. Report parsers → CanonicalFinding. |
| `compare` | 6 | VCF-level: `bcftools norm` + hap.py / `rtg vcfeval` with stratification (`--happy-gender` passes the sample's sex as hap.py `--gender`, needed for haploid chrX/chrY truth). Finding-level: matching expected vs observed findings; failure-stage attribution from intermediates (true alignments, BAM coverage, raw vs filtered VCF, annotated VCF, report). |
| `metrics` | – | Per category / variant type / stratum; Wilson intervals; acceptance criteria from config. |
| `diff` | – | Regression diff of two runs: new FN/FP, classification changes, stage changes. |
| `fidelity` | – | Simulated vs real data (HG002): planting summary across run directories; `samtools stats` / `bcftools stats` in the QC regions; hap.py (vcfeval) with stratifications; recall gaps against the real baseline. Steps write raw tool outputs; `fidelity.json` is assembled from whatever is present. |
| `report` | – | JSON for machines, HTML/Markdown for humans, figures used in the paper. |
| `schemas` | – | Pydantic models implementing `schemas.md`. |
| `genome.refdiff` | – | `clinval check-references`: the suite (simulation) reference and the pipeline's reference must share the primary contigs; BED of bases masked only in the simulation reference (GIAB-masked GRCh38, so reads of a GRCh38 false duplication come only from the true copy). |
| `genome` | – | Reference access and contig-name resolution (`1` → `chr1`), BED interval sets, variant-type classification, left normalisation, genomic HGVS (`NC_…:g.`) → VCF. |
| `config` | – | Suite YAML loading (duplicate keys rejected), path resolution relative to the suite file, input checks. |

## Variant selection (step 1)

1. Each source is read once; every record is evaluated against every rule of that source.
   ClinVar and GIAB VCFs are streamed with pysam. The ClinGen eRepo export (≈10^4 rows) is read
   into memory, mapped to GRCh38 from its genomic HGVS (checked against the FASTA and
   left-normalised) and merged per allele before it is offered to the rules (`docs/schemas.md`,
   *ClinGen eRepo*). Records with non-simple alleles (multi-allelic, symbolic, `.`, `N`) or contigs absent
   from the reference are skipped and counted as warnings.
2. Each rule applies a chain of filters and counts the candidates remaining after each one
   (`selection_report.md`). Passing candidates go to a *bucket* (gene, ladder bin, variant type,
   requested id or `total`) that keeps a bounded pool of the best-ranked candidates
   (requested + max(1000, 10 × requested)); an exhausted truncated pool is an error, never a
   silent shortfall.
3. **Ranking**: ClinVar review stars descending, then `derive_int(seed, "rank", variant_id)`
   ascending; eRepo and GIAB: hash only (all eRepo records are expert-panel classifications);
   explicit variants: listed order.
4. **Deduplication** is suite-wide: rules in YAML order, buckets in configuration order,
   candidates in rank order. A candidate is rejected (reason logged) if its REF span overlaps an
   accepted variant or the gap between REF spans is below `min_distance_bp`; the next candidate
   is tried. Variants in different samples are also kept apart.
5. Shortfalls are reported per bucket. `required` rules that are unsatisfied and samples
   without variants make the command exit with status 2 after writing all outputs.

## Truth encoding (step 2)

REF check against the FASTA for all variants (all mismatches reported at once, exit 3) →
genotype assignment per sample (`truth/genotype.py`) → optional background merge
(`truth/background.py`: background VCF streamed per contig, records near planted variants
removed and logged, unphased hets phased deterministically) → pysam VCF writing (planted and
background records merged by position) → `bcftools norm` (must be a no-op for planted records,
exit 3 otherwise; background multi-allelic records are split) → bgzip + tabix (pysam) →
expected outcomes and sample spec.

## Read simulation (step 3)

`prepare-sim` / `simulate` read the run directory written by `build-truth` (the truth VCF must
match the SHA-256 recorded in the manifest). Mason's input constraints were established from its
source (`apps/mason2`, seqan-v2.4.0) and test runs and are
enforced by `simprep` instead of being assumed:

1. Segments per sample (`docs/schemas.md`, *Simulation outputs*). Mason takes the number of
   haplotypes from the first VCF record and uses it for every contig, so diploid and haploid
   segments (male chrX/chrY outside the PARs) are simulated in **separate runs**, with fragments
   scaled by ploidy (haploid loci at half depth).
2. Planted variants are reserved on their haplotypes first; conflicts among them are errors
   (exit 4). The truth VCF is then streamed: records are shifted to segment coordinates,
   MNV/delins decomposed (Mason would silently skip them), GT-only records written, background
   records that cannot be represented dropped and logged.
3. `mason_simulator` per group, or per batch of it with `simulation.mason_batch_bp` (runner
   `mason_simulator`, cwd = `<sample>.sim/`, `TMPDIR=<sample>.sim/` because Mason keeps all reads in
   temporary files until its final join), read names prefixed `<sample>.<run>.`; seed derived from
   the suite seed and the run name; `threads`/`chunk_size` from the profile because they change
   the output. Each run is merged into the outputs (in a background thread, while Mason simulates
   the next run) and its raw files deleted, so the disk peak stays near that of the largest run.
   The lift edits the SAM text of each record (RNAME, POS, RNEXT, PNEXT) and parses it back
   against the genome header (htslib), ~25 % faster than via `to_dict`/`from_dict`.
4. FASTQ concatenated per mate into one gzip stream (ISA-L level 3 via `isal.igzip_threaded`,
   1 MiB blocks, no file name, `mtime=0`); true alignments lifted to genome coordinates (contig,
   position, mate) and sorted with `samtools sort --no-PG`. Compression, the lift's BGZF I/O,
   `samtools sort` and `samtools index` use the profile's `threads`; their outputs do not depend
   on the thread count. The gzip bytes depend on the ISA-L level and the `isal` version.
5. Planting check: every planted variant must be carried by ≥ `min_alt_reads` reads on each
   carrying haplotype (`oH` tag); VAF and per-haplotype counts are kept for the comparison (`not_simulated`
   / `not_covered`). Incomplete ⇒ exit 4 after writing all outputs.

Region mode simulates only planted loci (fast; reads from paralogues come only from
`extra_regions_bed`); genome mode simulates whole chromosomes (realistic off-target mapping, ~30x
WGS cost). Mason places no fragments in N gaps but distributes fragments by total length, so
contigs with long gaps get slightly higher depth in their non-N part (documented, not corrected).

## CLI and exit codes

`clinval validate-suite`, `select`, `build-truth`, `prepare-sim`, `simulate`, `validate-pipeline`, `run-sut`, `fidelity`, `export-schemas`.
`build-truth` and `simulate` take `--runtime` (RuntimeConfig). Exit codes: 0 success;
1 invalid suite, missing input or tool failure; 2 selection unsatisfied; 3 truth encoding failed
(REF mismatch, including reference bases stated in eRepo HGVS); 4 simulation failed (planted
variant cannot be represented for Mason, or planting check incomplete).
Logs go to `<run dir>/clinval.log`.

## Reference data

`scripts/fetch_data.py` (datasets in `scripts/data_sources.yaml`) is the only networked step:
download to `*.part`, SHA-256 verified against the pinned value, then renamed; unpinned files
(moving exports such as the ClinGen eRepo) only with `--allow-unpinned`; preparation steps
(reference decompression + `.fai`); report `data/fetch_report.json`. Datasets without a public URL
(eRepo v2.2.0) are listed with instructions. `--ca-bundle` (or `SSL_CERT_FILE`) for
TLS-inspecting proxies.

## Cross-cutting

- **Subprocess wrapper** (`insilico_clinval/_exec.py`): no shell, return-code checks, stderr
  capture, wall time and peak RSS of the child (`os.wait4`). Pipeline tools are called only
  through `Toolbox` (`insilico_clinval/_runner.py`), which builds the runner command, captures the
  tool version (`<tool> --version`, first line, through the same runner) and records version,
  container digest and step in the manifest.
- **Provenance** (`insilico_clinval/provenance.py`): run manifest create/extend/write, streamed
  SHA-256, framework version, git commit and dirty flag.
- **Determinism** (`insilico_clinval/rng.py`): every random decision is
  `int.from_bytes(sha256(f"{seed}|{purpose}|{key}")[:8], "big")`, independent of iteration
  order. Purposes: `rank` (key = variant id), `haplotype` (key = `{sample}|{variant_id}`),
  `phase_set` (key = `{sample}|{chrom}:{PS}`), `bg_haplotype` (key = `{sample}|{chrom}:{pos}`,
  unphased background hets). JSON is written with sorted keys, all lists are
  sorted, `bcftools` runs with `--no-version`, output paths are relative to the run directory;
  timestamps and timings exist only in `run_manifest.json` and `clinval.log`. Byte identity holds for
  the same tool versions (recorded in the manifest): e.g. the true-alignment BAM written by
  samtools 1.16.1 (container) and 1.23 (local) contains identical records but different BGZF bytes.
  Mason's output also depends on `threads` and `chunk_size` (profile).
- **Constants** (`insilico_clinval/constants.py`): GRCh38 PAR coordinates (GRC); RefSeq
  accessions of GRCh38 chromosomes (NCBI assembly report GRCh38.p14).
- **Orchestration**: every step is a `clinval` sub-command, run in sequence and resumable step by
  step (the run manifest binds a run directory to its suite); `examples/quickstart/run_quickstart.sh`
  chains them. There is no all-in-one `run-suite` command yet. **No Docker Compose.**
- **Execution / runner abstraction** (`insilico_clinval/_runner.py`): every external tool is invoked
  through a runner selected per tool in a RuntimeConfig (`--runtime`, `configs/runtime/`; separate
  from the suite so one suite runs locally or in containers): `docker` (`docker run --rm --user
  UID:GID` of an image pinned as `name@sha256:…`; the directories a call needs are bind-mounted at
  identical resolved paths, nested mounts collapsed, working directory set) or `local` (tool on
  PATH, the default). For docker steps the manifest records the image digest and no peak RSS (the
  measured process is the docker client). An `apptainer` runner may be added later for HPC. Tool images live in `containers/<tool>/Dockerfile`:
  `clinval` (framework + bcftools/samtools/tabix), `mason` (built from SeqAn `seqan-v2.4.0`;
  Debian's package lacks `mason_simulator`), `happy` (hap.py + rtg vcfeval, bioconda, amd64),
  `bias2015` (built from the upstream repository at a pinned tag; annotation data in a mounted
  directory, not in the image), `nirvana` (Nirvana 3.18.1). BIAS-2015's VEP mode uses the
  upstream `ensemblorg/ensembl-vep:release_115.0` image unchanged.
- **Nextflow / nf-core/sarek runs on the host**, not in a container: `nextflow run nf-core/sarek
  -r <pinned> -profile docker`; Sarek pulls its own per-process containers. The secondary adapter
  calls the host `nextflow` binary (runner `local`) — this avoids docker-in-docker and path issues.
