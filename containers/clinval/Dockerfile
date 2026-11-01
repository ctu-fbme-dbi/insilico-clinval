# Framework image: clinval + bcftools/samtools/tabix (Debian bookworm packages, pinned versions).
# Build from the repository root:  docker build -f containers/clinval/Dockerfile -t clinval:dev .
FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258
RUN apt-get update && apt-get install -y --no-install-recommends \
        bcftools=1.16-1 samtools=1.16.1-1 tabix=1.16+ds-3 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/insilico-clinval
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .
LABEL org.opencontainers.image.title="clinval"
ENTRYPOINT []
