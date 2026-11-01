# Running the ACMG SF v3.2 test sets on a Linux server

Step-by-step manual for producing the 10 paired-end FASTQ sets (5 suites × `positive`/`negative`)
on a Linux server. Steps 1–6 are one-time setup; steps 7–9 run and check the simulations. Running
the pipelines under test and scoring them is described in [reproducing.md](reproducing.md).

## The test sets

Five suites in `configs/suites/acmg_sf32_*.yaml` plant ClinGen eRepo (VCEP-curated) variants in
ACMG SF v3.2 genes (Miller et al., Genet Med 2023;25(8):100866), restricted to SF-compatible
diseases (MONDO allowlist), into the GIAB HG002 background. Mason simulates 2×150 bp reads at 30x
in genome mode.

Each suite has two samples:

- **`positive`**: findings that must be reported as secondary findings (P/LP in AD/SD/XL genes,
  and biallelic P/LP in AR genes).
- **`negative`**: findings that must not be reported (B/LB, VUS, and single heterozygous P/LP in
  AR genes).

| Suite | Simulated contigs | Planted variants |
|---|---|---|
| `acmg_sf32_short1` | chr17, chr19 | 33 |
| `acmg_sf32_short2` | chr12–chr15 | 30 |
| `acmg_sf32_short3` | chr1, chr10, chrX | 23 |
| `acmg_sf32_full1` | chr1–22, chrX, chrY | 88 (seeded random draw) |
| `acmg_sf32_full2` | chr1–22, chrX, chrY | 69 (one P/LP per SF gene) |

The article's `acmg_sf32_full1` run (simulated 2026-10-22, seed 20261018) was built from
`configs/suites/acmg_sf32_full1.prev.yaml`, the suite as of commit `baa1de2` (sha256
`81c39442…`). Commit `d0e9f96` later narrowed each rule's MONDO allowlist in
`acmg_sf32_full1.yaml`; the author kept the simulated run. Every later command on that run
directory (`run-sut`, comparison) must use `acmg_sf32_full1.prev.yaml`: the run manifest accepts
only the config it was written for. Do not edit the `.prev` file, not even its comments.

## 1. Check the server has the tools

```bash
python3 --version          # needs >= 3.10
bcftools --version
samtools --version         # tabix comes with htslib
docker info                # must work without sudo (your user in the "docker" group)
df -h /path/to/big/volume  # plan on ~0.8–1 TB free for all 10 pairs
```

If bcftools or samtools are missing, install them with either of:

```bash
sudo apt install bcftools samtools tabix
conda install -c bioconda bcftools samtools
```

Disk estimate:

- The short sets need tens of GB.
- A 30x genome is roughly 70–100 GB of FASTQ.gz, plus a ~75 GB `true.bam`, plus a few hundred GB
  of scratch space while Mason runs.

## 2. Clone the repo and install

```bash
git clone <repository-url> clinval
cd clinval
python3 -m venv .venv
source .venv/bin/activate          # repeat in every new shell
pip install -e ".[dev]"
clinval --help
```

## 3. Build the Mason image

This takes about a minute and compiles Mason from the SeqAn 2.4.0 sources:

```bash
docker build -t clinval-mason:2.4.0 containers/mason
```

## 4. Send run outputs to the large disk

The suites write to `runs/acmg_sf32_*`. `runs/` is git-ignored, so make it a symlink to the big
volume:

```bash
mkdir -p /mnt/data/projects/clinval-runs
ln -s /mnt/data/projects/clinval-runs runs
```

## 5. Check the installation

This takes seconds and needs no downloads.

```bash
pytest -m "not integration" -q

clinval build-truth --suite tests/fixtures/suite_erepo.yaml --out runs/fixture_erepo
clinval simulate    --suite tests/fixtures/suite_erepo.yaml --out runs/fixture_erepo \
                    --runtime configs/runtime/docker-local-build.yaml

cat runs/fixture_erepo/positive.planting_check.tsv
# every line should end with "ok"
```

If this works, Docker, Mason, bcftools and samtools are all wired correctly.

## 6. Get the input data

The reference genome and the HG002 background have fixed checksums. They are about a 1 GB
download, and 3.1 GB once the reference is unpacked.

```bash
python scripts/fetch_data.py --only grch38,giab_hg002_v4.2.1
```

The eRepo export must be **this exact file**, not a fresh download. eRepo changes daily, and the
suites' variant selection was tuned on the 2026-10-15 export (archived with the article's data).
Copy it to `data/clingen/erepo_2026-10-15.tsv` and check it:

```bash
echo "a7aba27b6308da9f3c8b64374c236abc6681f1bcb25d448fb22009e324899795  data/clingen/erepo_2026-10-15.tsv" | sha256sum -c -
# expected: data/clingen/erepo_2026-10-15.tsv: OK
```

## 7. Validate the suites

```bash
for s in short1 short2 short3 full1 full2; do
  clinval validate-suite configs/suites/acmg_sf32_$s.yaml
done
```

This checks that every input file exists, so a missing download or wrong path shows up now rather
than hours into a run.

## 8. Run the simulations

### Use tmux

Run inside `tmux` so the job survives an SSH disconnect:

| Action | Command |
|---|---|
| Start a session | `tmux new -s sim` |
| Detach | `Ctrl-b d` |
| Reattach | `tmux attach -t sim` |

### Short sets first

Each short set covers a few chromosomes and takes minutes to an hour.

```bash
source .venv/bin/activate
RT=configs/runtime/docker-local-build.yaml
for s in short1 short2 short3; do
  clinval build-truth --suite configs/suites/acmg_sf32_$s.yaml || break
  clinval simulate    --suite configs/suites/acmg_sf32_$s.yaml --runtime $RT || break
done
```

### Then the full sets

Run one sample at a time to limit peak disk use:

```bash
for s in full1 full2; do
  clinval build-truth --suite configs/suites/acmg_sf32_$s.yaml || break
  clinval simulate    --suite configs/suites/acmg_sf32_$s.yaml --runtime $RT --sample positive || break
  clinval simulate    --suite configs/suites/acmg_sf32_$s.yaml --runtime $RT --sample negative || break
done
```

`|| break` stops at the first failure instead of carrying on.

### Logging and progress

`clinval` prints only a one-line summary (stdout) or an `error: ...` message (stderr) per command;
the detail goes to files in the run directory:

| File | Written | Contains |
|---|---|---|
| `clinval.log` | appended by every command, while it runs | selection decisions (rejected variants, skipped eRepo alleles) and, for `simulate`, progress of every stage with elapsed times |
| `run_manifest.json` | at the end of each command (also on a tool error) | every step with wall time, peak memory and exit code; tool versions; input/output SHA-256; warnings |
| `<sample>.sim/mason_<group>.log` | after Mason finishes | Mason's own output |
| `<sample>.planting_check.tsv` | at the end of `simulate` | per planted variant: reads, VAF, `ok` or the failure reason |

`select`, `build-truth`, `prepare-sim` and `simulate` print their phases and progress on the
console (stderr, prefixed with `HH:MM:SS`); the full log, including per-candidate selection
details, is in `clinval.log`. `build-truth` prints per sample:

```text
selecting variants: <n> rule(s)
selected <n> variants in <H:MM:SS>
truth: <n> selected variants, 2 sample(s)
positive: writing truth VCF: <n> planted variants + GIAB_HG002 background
positive: truth VCF: 1,000,000 records written (at chr4)     # every million records
positive: truth VCF written in <H:MM:SS>: <n> records
positive: normalising (bcftools norm), compressing and indexing
positive: truth VCF normalised and indexed in <H:MM:SS>
positive: truth finished in <H:MM:SS>
```

On the full genome (≈ 4 M HG002 background records) this takes about a minute per sample.
Follow a running `simulate` from a second shell (or after reattaching) with:

```bash
tail -f runs/acmg_sf32_<set>/clinval.log
```

Per sample it logs these lines (timestamps omitted; `<...>` are the actual values):

```text
positive: preparing Mason inputs (genome mode)
positive: Mason inputs ready in <H:MM:SS>: <n> segments, <n> records (diploid <n>, haploid <n>), <n> background dropped
positive: simulating <k> Mason run(s), 30x, 2 x 150 bp, 8 Mason thread(s), <n> alignments
positive: mason_simulator diploid-1 started: <n> fragments over <n> bp   # full1/full2: batches
positive: Mason output: runs/acmg_sf32_<set>/positive.sim/mason_diploid-1.log
positive: mason_simulator diploid-1 done in <H:MM:SS>
positive: mason_simulator diploid-2 started: ...         # next Mason run starts at once ...
positive: diploid-1: FASTQ appended in <H:MM:SS>         # ... while run 1 is merged alongside
positive: lifted 10,000,000 / <n> alignments (<p>%)      # every 10 million, running total
positive: diploid-1: <n> true alignments lifted in <H:MM:SS> (<so far> / <n>)
...                                                      # the same for every run
positive: mason_simulator haploid-1 started: ...         # chrX/chrY outside the PARs
positive: waiting for the last merge (FASTQ + lift)
positive: sorting and indexing the true BAM (samtools)
positive: true BAM sorted and indexed in <H:MM:SS>
positive: planting check <ok>/<n> variants ok
positive: simulation finished in <H:MM:SS>
```

Without `simulation.mason_batch_bp` (the short suites) there is one run per ploidy group,
named `diploid` and `haploid`.

Mason's own output (stderr) is written live to `<sample>.sim/mason_<group>.log` (the path is
logged as `<sample>: Mason output: …`), not to the console, so the "started" line is the last
console line for the length of the Mason run. Follow Mason's progress (one dot per chunk, per
contig and haplotype) with `tail -f runs/acmg_sf32_<set>/<sample>.sim/mason_diploid.log`. The log
is kept when Mason fails, e.g. with a segmentation fault (exit status 139). All Mason inputs are
prepared first, for every sample of the command, before the first Mason run starts.

To also keep the console output with timestamps (bash), use this instead of the plain loop above:

```bash
mkdir -p logs
for s in short1 short2 short3; do
  ( set -e
    date; clinval build-truth --suite configs/suites/acmg_sf32_$s.yaml
    date; clinval simulate    --suite configs/suites/acmg_sf32_$s.yaml --runtime $RT
    date ) 2>&1 | tee logs/$s.console.log
  [ "${PIPESTATUS[0]}" -eq 0 ] || break   # exit status of the subshell, not of tee
done
```

`logs/` is git-ignored.

During the silent Mason stage:

```bash
docker ps --filter ancestor=clinval-mason:2.4.0   # Mason still running?
tail -f runs/acmg_sf32_<set>/<sample>.sim/mason_diploid.log   # contig being simulated, dots
ls -lh runs/acmg_sf32_<set>/<sample>.sim/          # Mason's FASTQ/BAM files growing
```

### Notes

- **Disk for Mason's temporary files.** Mason keeps all simulated reads in temporary files
  (created in `$TMPDIR`, unlinked at once, so invisible to `ls`) until it joins them into its
  outputs at the end. `simulate` sets `TMPDIR` to `<sample>.sim/`, so they land on the run disk.
  At the join, temporary files and outputs coexist. Measured on chr4 (19 M fragments): 1.45 KB of
  temporary files and 0.97 KB of Mason output per fragment, peak 2.4 KB per fragment; a full 30x
  genome in one Mason run (288 M fragments) therefore peaks at about 700 GB on the run disk. The
  full suites set `simulation.mason_batch_bp: 500000000` (8 Mason runs); each run's raw outputs
  are deleted as soon as it is merged (while the next run is simulated, so one run's raw outputs,
  ~45 GB, may coexist with the next run). Measured on a full `acmg_sf32_full1` sample (2026-10-17,
  Apple M4, 8 threads): peak **317 GB** on the run disk, during the final `samtools sort` (the
  largest batch, chr13–chr17, peaked at 259 GB including the outputs merged so far); final
  outputs 143 GB per sample (2 × 38.5 GB FASTQ.gz, 66 GB true BAM); peak memory 7.5 GB
  (`samtools sort`); `simulate` 3:25 h (Mason 2:21, lift 0:44, sort + index 0:15). Plan on
  ≥ 350 GB free for the sample being simulated plus 145 GB for every finished one. If the disk
  fills, Mason dies with a segmentation fault (exit status 139) because it does not check write
  errors.

- **Run time.** A full 30x sample took 3:25 h locally (Apple M4, 8 threads; Mason 2:21, lift 0:44,
  sort + index 0:15); on a Linux server Mason was faster and the lift slower. Mason is the slowest part. FASTQ compression (ISA-L), the lift of Mason's
  alignments back to genome coordinates (`simulate/merge.py`) and `samtools sort` use the
  profile's threads; the lift loop itself is still one Python thread.
- **Parallel runs.** Different sets can run at the same time in separate tmux windows, because each
  writes to its own directory. That needs enough cores (8 threads each) and disk. Never run two
  commands on the *same* set at once.
- **Do not change the thread count.** Reads are reproducible only for the same seed and thread
  count; the profile `illumina_2x150_30x_t8.yaml` uses 8.
- **SLURM.** Wrap each `build-truth` + `simulate` pair in an `sbatch` script requesting 8 CPUs.
  Docker must be available on the compute nodes.

## 9. Check the results

### Exit codes

| Code | Meaning |
|---|---|
| 0 | OK |
| 1 | Invalid input or tool error |
| 2 | Selection not satisfied |
| 3 | Truth encoding failed |
| 4 | Simulation failed, or a planted variant has no reads |

### Commands

```bash
# planting complete for every sample (all should say true)
grep -H '"planting_complete"' runs/acmg_sf32_*/*.simulation.json

# any planted variant not confirmed in the reads (should print nothing)
awk -F'\t' 'FNR>1 && $NF!="ok" {print FILENAME": "$0}' runs/acmg_sf32_*/*.planting_check.tsv

# simulated contigs, e.g. only chr17 and chr19 for short1
cat runs/acmg_sf32_short1/positive.sim/segments.tsv

# the 10 FASTQ pairs
ls -lh runs/acmg_sf32_*/*.R[12].fastq.gz
```

### Deliverables

Each `runs/acmg_sf32_<set>/` contains these files, for both `positive` and `negative`:

| File | Purpose |
|---|---|
| `<sample>.R1.fastq.gz`, `<sample>.R2.fastq.gz` | Input for the pipeline under test |
| `<sample>.truth.vcf.gz` | What should be called |
| `<sample>.expected_outcomes.json` | What should be classified and reported |
| `<sample>.true.bam` | Mason's true read alignments, for debugging |

## Troubleshooting

| Symptom | Fix |
|---|---|
| `validate-suite` reports a missing file | Repeat step 6 and check the paths under `data/`. |
| "run_manifest.json was written for a different suite config or seed" | The suite was edited after a run. Delete `runs/acmg_sf32_<set>` and rerun that set. |
| Exit 4, some variants `no_reads` | Inspect `<sample>.planting_check.tsv`. This should not happen at 30x. |
| Disk full during a full set | Delete the failed sample's `<sample>.sim/`, free space, and rerun only that `--sample`. |

## Known limitations

- **Gene coverage.** Only 29 of the 81 SF v3.2 genes have eRepo curations. Not covered: TTN, LMNA,
  SCN5A, KCNH2, APOB, PCSK9, CALM1–3, HFE, MUTYH, ATP7B, GLA and others.
- **Background in short sets.** The truth VCF still holds HG002 background for all of chr1–22.
  Subset it to the simulated contigs before comparing calls.
- **chrX/chrY.** HG002 v4.2.1 covers chr1–22 only, so chrX/chrY carry planted variants only.
- **Held-out flag.** The suites configure no eRepo v2.2.0 reference, so `held_out` in
  `variants.selected.tsv` is false for every variant. `clinval held-out` computes the real flag
  afterwards (`held_out.tsv`, see [reproducing.md](reproducing.md)).

## Next step

Run the pipelines under test on the FASTQs and score them: [reproducing.md](reproducing.md),
section 4.
