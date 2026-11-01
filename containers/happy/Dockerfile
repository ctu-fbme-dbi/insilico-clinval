# hap.py + rtg-tools (VCF comparison) from bioconda. Bioconda builds hap.py for
# linux-64 only; on arm64 hosts the image runs under emulation (--platform linux/amd64).
#   docker build --platform linux/amd64 -t clinval-happy:0.3.15 containers/happy
FROM --platform=linux/amd64 mambaorg/micromamba:2.0.5@sha256:8d0d4cad074afe3ac93b58dce88bd2659153d3fa6de3a6968f501710e8d3c832
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        hap.py=0.3.15 rtg-tools=3.12.1 \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
LABEL org.opencontainers.image.title="happy" org.opencontainers.image.version="0.3.15"
ENTRYPOINT []
