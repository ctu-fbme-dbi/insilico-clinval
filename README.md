# insilico-clinval

End-to-end in silico validation and regression testing of clinical genome sequencing pipelines,
from simulated reads to reportable findings.

Code companion of the article *"An end-to-end in silico framework for validation and regression
testing of clinical genome sequencing pipelines, from simulated reads to reportable findings"*
(submitted to Frontiers in Genetics). `insilico-clinval` is a working name; the command is
`clinval`.

## What it does

You give it a **test suite** (which variants to plant, in which genome, at what depth) and a
**pipeline config** (how to run the pipeline you want to test). It then:

```
suite.yaml ─► 1 select ─► 2 build-truth ─► 3 simulate ─► 4-5 run-sut ─────────► 6 compare ─► metrics / diff / report
              variants    truth VCF +       Mason reads   your pipeline:          outcome per planted variant,
                          expected          (FASTQ) +     FASTQ ─► VCF ─► ACMG    failure stage, false reports
                          outcomes          true BAM      report (black box)
```

1. **Select** test variants: expert-curated ClinGen eRepo variants (pathogenic/likely pathogenic
   and, as negative controls, benign/likely benign in the same genes), ClinVar, GIAB, or your own
   synthetic variants.
2. **Build the truth**: a phased truth VCF on the GIAB HG002 background, plus the *expected
   outcome* of every planted variant (its expected ACMG class and whether it must be reported).
3. **Simulate** paired-end reads de novo from the modified diploid genome with Mason (SeqAn), and
   check that every planted variant is really in the reads.
4. **Run the system under test**: secondary analysis (FASTQ → VCF, e.g. nf-core/sarek) and
   tertiary analysis (VCF → classified, reported findings, e.g. BIAS-2015), as black boxes.
5. **Compare** with the truth: one outcome per planted variant (correct, not called,
   misclassified, not reported, ...), the pipeline **stage** where each failure happened, and
   every false report.
6. **Score and diff**: sensitivity and specificity of reporting with 95 % Wilson intervals,
   acceptance criteria, and a variant-by-variant regression diff between two pipeline versions.

The framework is meant for **supplemental validation and change-control (regression) testing**,
not as a replacement for validation with human samples (Roy et al., 2018, Recommendations 14
and 17).

## Quickstart (about 10 minutes, no downloads)

Runs all six steps on a tiny synthetic genome bundled in `tests/fixtures/`. Works on Linux and
macOS (Intel or Apple Silicon).

### Requirements

| Tool | Check | Install |
|---|---|---|
| Python ≥ 3.10 | `python3 --version` | |
| bcftools, samtools | `bcftools --version` | `brew install bcftools samtools` · `sudo apt install bcftools samtools tabix` · `conda install -c bioconda bcftools samtools` |
| Docker | `docker info` | Docker Desktop or Docker Engine (your user must be able to run `docker` without `sudo`) |

### 1. Install

```bash
git clone <repository-url> insilico-clinval
cd insilico-clinval
python3 -m venv .venv
source .venv/bin/activate          # in every new shell
pip install -e ".[dev]"
clinval --help
```

(With [uv](https://docs.astral.sh/uv/): `uv venv && uv pip install -e ".[dev]"`.)

### 2. Build the two tool images

```bash
docker build -t clinval-mason:2.4.0 containers/mason                          # ~1 min
docker build --platform linux/amd64 -t clinval-happy:0.3.15 containers/happy  # ~5 min
```

Mason is compiled from the SeqAn 2.4.0 sources; hap.py comes from bioconda (amd64 only; on Apple
Silicon Docker emulates it).

### 3. Run the end-to-end example

```bash
examples/quickstart/run_quickstart.sh          # writes runs/quickstart/, takes ~1 min
```

The script prints every `clinval` command it runs, so you can repeat or change any step by hand.
It tests a deliberately simple pipeline, `toy` (bcftools calls variants on Mason's true
alignments; a lookup-table classifier labels known variants and calls everything else a VUS), and
a modified version of it, `toy-noPP5BP6`, whose classifier ignores its knowledge base.

### What to look at

| File in `runs/quickstart/` | What it shows |
|---|---|
| `report.md` | metrics with 95 % CIs, acceptance criteria, and the regression diff of `toy-noPP5BP6` against `toy` |
| `sut/toy/compare/outcomes.tsv` | one row per planted variant: `outcome_class`, the failure `stage`, expected vs observed class |
| `selection_report.md` | which variants were selected and why |
| `positive.planting_check.tsv` | every planted variant confirmed in the simulated reads |
| `run_manifest.json` | provenance: versions, checksums of all inputs and outputs, wall time per step |

The toy pipeline fails on purpose in two places, and the framework finds both:

- `curated:chr2:600:GAG:TT` → `not_called`, stage `mapped`: the reads were aligned correctly, but
  the caller split the delins into two records (a calling failure);
- `synthetic_novel:chr2:5100:G:T` → `misclassified`, stage `annotated`: the variant was called,
  but the classifier knows only previously classified variants, so this novel loss-of-function
  variant became a VUS and was not reported (an interpretation failure).

The regression diff marks the candidate `toy-noPP5BP6` as a **regression** (11 new failures), and
the acceptance criteria of `configs/acceptance/article-default.yaml` fail for planted-variant
calling and report sensitivity.

### Run the tests

```bash
pytest -m "not integration"                                    # ~250 unit tests, seconds
CLINVAL_TEST_MASON_IMAGE=clinval-mason:2.4.0 \
CLINVAL_TEST_HAPPY_IMAGE=clinval-happy:0.3.15 pytest           # + integration tests (bcftools, samtools, docker)
ruff check . && ruff format --check . && mypy                  # code quality
```

An opt-in test runs nf-core/sarek on the fixture genome (`CLINVAL_TEST_SAREK=1`, needs Nextflow).

## Test your own pipeline

A pipeline is described by a YAML **pipeline config** with a `secondary` and a `tertiary`
section. Each section uses an adapter:

| Section | Adapters |
|---|---|
| `secondary` (FASTQ → VCF) | `sarek` (nf-core/sarek, run by the framework or adopted with `--sarek-results`), `command` (any command line or container), `precomputed` (VCFs produced elsewhere, e.g. DRAGEN) |
| `tertiary` (VCF → findings) | `bias2015` (BIAS-2015 with Nirvana or VEP), `command` (any classifier whose TSV/CSV/JSON output is mapped by a parser), `none` |

Start from `examples/quickstart/toy_pipeline.yaml.in` (command adapters, the simplest case),
`configs/adapters/example_blackbox_classifier.yaml` (Sarek + a closed classifier), or the
article's configs in `configs/adapters/`. Then:

```bash
clinval validate-pipeline my_pipeline.yaml
clinval run-sut --suite SUITE.yaml --out runs/X --pipeline my_pipeline.yaml
clinval compare --suite SUITE.yaml --out runs/X --pipeline my_pipeline.yaml --runtime configs/runtime/docker.yaml
```

Details: [User manual](docs/manual.md), sections 7–8.

## Real data and the article

| Goal | Read |
|---|---|
| Run on GRCh38 with real ClinGen eRepo variants (a small demo suite, minutes) | [User manual](docs/manual.md), sections 5–6 |
| Reproduce the article's test sets and results | [Reproducing the article](docs/reproducing.md) |
| Run the 30x whole-genome suites on a Linux server | [ACMG SF v3.2 test sets on a server](docs/acmg_sf32_server.md) |

Reference data are downloaded by `scripts/fetch_data.py` and verified by SHA-256; nothing large
is stored in the repository. The article's runs need roughly 1 TB of disk and a Linux server
with Docker.

## Documentation

| Document | Content |
|---|---|
| [docs/manual.md](docs/manual.md) | user manual: installation, data, every command, outputs, configuration, exit codes, troubleshooting |
| [docs/reproducing.md](docs/reproducing.md) | which suite, config and command produced each result of the article |
| [docs/acmg_sf32_server.md](docs/acmg_sf32_server.md) | running the large suites on a server (disk, time, monitoring) |
| [docs/architecture.md](docs/architecture.md) | design: data flow, modules, determinism, provenance |
| [docs/schemas.md](docs/schemas.md), [docs/json_schema/](docs/json_schema/) | data contracts of every input and output file |
| [containers/README.md](containers/README.md) | tool images |

## Repository layout

```
src/insilico_clinval/  the framework (CLI: cli.py)
  variants/ truth/ simprep/ simulate/   steps 1-3: selection, truth, read simulation
  adapters/                             steps 4-5: secondary and tertiary analysis
  compare/ metrics/ diff/ report/       step 6: comparison, metrics, regression diff
  fidelity/                             simulated vs real genome (article, Results 3.1)
  schemas/                              pydantic models of all data contracts
configs/
  suites/      test suites (the article's ACMG SF v3.2 sets, difficult contexts, supplementary experiments, demo)
  adapters/    pipeline configs (Sarek + BIAS-2015 variants used in the article, templates)
  profiles/    read-simulation profiles      runtime/  where external tools run (local or docker)
  acceptance/  acceptance criteria           regions/  BED files     fidelity/  fidelity experiment
  synthetic/   inputs of scripts/generate_synthetic.py
containers/    Dockerfiles: Mason, hap.py, Nirvana, BIAS-2015, the framework itself
examples/      quickstart: toy pipeline and end-to-end script
scripts/       fetch_data.py (checksummed downloads), generate_synthetic.py, background_diff.py
tests/         unit/ (fixtures only), integration/ (external tools), fixtures/ (tiny synthetic genome)
```

## Citation and licence

Please cite the article (reference to be added on publication) and the archived release;
`CITATION.cff` holds the software citation. Release 1.0.0 is the version described in the article and
is archived on Zenodo: https://doi.org/10.5281/zenodo.23859473 (all versions:
https://doi.org/10.5281/zenodo.23859472). Changes are listed in `CHANGELOG.md`.

The code is source-available under the Academic and Research License 1.0 (`LICENSE.md`): academic
research, teaching and noncommercial research are free. Clinical, diagnostic and
laboratory-operations use (including validation of clinical pipelines, also in public hospitals) and
commercial use require a separate licence on request; see `COMMERCIAL-LICENSE.md`. External tools
and data sets keep their own licences; BIAS-2015 (AGPL-3.0) is called as an external tool and is
not part of this repository.
