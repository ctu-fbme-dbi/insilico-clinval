# BIAS-2015 (Bitscopic, AGPL-3.0 for noncommercial academic use; commercial use needs a Bitscopic
# licence) cloned from the upstream repository at a pinned tag during the build. Its code is
# never copied into this repository. Data (hg38 bundle) is mounted at run time.
#   docker build -t clinval-bias2015:v3.0.0 containers/bias2015
FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258
ARG BIAS_TAG=v3.0.0
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates \
    && git clone --depth 1 --branch ${BIAS_TAG} https://github.com/bitscopic/BIAS-2015.git /opt/BIAS-2015 \
    && git -C /opt/BIAS-2015 rev-parse HEAD > /opt/BIAS-2015.commit \
    && rm -rf /opt/BIAS-2015/.git && apt-get purge -y git && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
# Wrappers: `bias_2015` (classifier; `--version` prints tag and commit) and `bias_2015_paths`
# (BIAS's own required_paths.json generator, which also points to data shipped in the repo).
RUN printf '#!/bin/sh\nif [ "$1" = "--version" ]; then echo "BIAS-2015 %s ($(cat /opt/BIAS-2015.commit))"; exit 0; fi\ncd /opt/BIAS-2015 && exec python3 bias_2015.py "$@"\n' "${BIAS_TAG}" \
       > /usr/local/bin/bias_2015 \
    && printf '#!/bin/sh\nif [ "$1" = "--version" ]; then echo "BIAS-2015 %s"; exit 0; fi\ncd /opt/BIAS-2015 && exec python3 src/scripts/create_new_required_paths_file.py "$@"\n' "${BIAS_TAG}" \
       > /usr/local/bin/bias_2015_paths \
    && chmod +x /usr/local/bin/bias_2015 /usr/local/bin/bias_2015_paths
LABEL org.opencontainers.image.title="bias2015" org.opencontainers.image.version="v3.0.0" \
      org.opencontainers.image.licenses="AGPL-3.0-only" \
      org.opencontainers.image.source="https://github.com/bitscopic/BIAS-2015/tree/v3.0.0"
ENTRYPOINT []
