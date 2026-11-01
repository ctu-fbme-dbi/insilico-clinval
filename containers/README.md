# Tool images

One image per external tool, selected per tool in a runtime config (`configs/runtime/`,
`clinval … --runtime FILE`; see `docs/architecture.md`, Execution). The image id (or digest) is
recorded in every run manifest. For published runs, push the images to a registry and pin them
by digest (`name@sha256:…`); a tag-only image produces a manifest warning.

| Directory | Image | Contents | Needed for | Build |
|---|---|---|---|---|
| `mason/` | `clinval-mason:2.4.0` | `mason_simulator` 2.0.9 from SeqAn `seqan-v2.4.0` (source build; Debian's `seqan-apps` lacks it) | `simulate` (quickstart) | `docker build -t clinval-mason:2.4.0 containers/mason` |
| `happy/` | `clinval-happy:0.3.15` | hap.py 0.3.15, rtg-tools 3.12.1 (bioconda, linux/amd64) | `compare`, `fidelity` (quickstart) | `docker build --platform linux/amd64 -t clinval-happy:0.3.15 containers/happy` |
| `nirvana/` | `clinval-nirvana:3.18.1` | Nirvana 3.18.1 on .NET 6; annotation data mounted at run time | BIAS-2015 with Nirvana | `docker build -t clinval-nirvana:3.18.1 containers/nirvana` |
| `bias2015/` | `clinval-bias2015:v3.0.0` | BIAS-2015 cloned from upstream at the pinned tag during the build (AGPL-3.0; never copied into this repository); data mounted at run time | BIAS-2015 | `docker build -t clinval-bias2015:v3.0.0 containers/bias2015` |
| `clinval/` | `clinval:dev` | the framework, bcftools 1.16, samtools 1.16.1, tabix | optional (`configs/runtime/docker.yaml`) | `docker build -f containers/clinval/Dockerfile -t clinval:dev .` |

BIAS-2015 in VEP mode uses the official `ensemblorg/ensembl-vep:release_115.0` image.
nf-core/sarek is not containerised here: Nextflow runs on the host with `-profile docker` and
pulls its own per-process containers.
