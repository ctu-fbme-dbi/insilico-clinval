# Mason (SeqAn 2.4.0) read simulator, built from the upstream release tag.
# Debian's seqan-apps package does not ship mason_simulator, hence the source build.
FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587 AS build
ARG SEQAN_TAG=seqan-v2.4.0
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates git cmake make g++ zlib1g-dev libbz2-dev python3 \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch ${SEQAN_TAG} https://github.com/seqan/seqan.git /src/seqan
WORKDIR /build
RUN cmake /src/seqan -DCMAKE_BUILD_TYPE=Release -DSEQAN_BUILD_SYSTEM=APP:mason2 \
        -DCMAKE_CXX_FLAGS="-w" \
    && make -j"$(nproc)" mason_simulator \
    && ./bin/mason_simulator --version

FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587
RUN apt-get update && apt-get install -y --no-install-recommends zlib1g libbz2-1.0 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /build/bin/mason_simulator /usr/local/bin/mason_simulator
LABEL org.opencontainers.image.title="mason" \
      org.opencontainers.image.version="2.4.0" \
      org.opencontainers.image.source="https://github.com/seqan/seqan/tree/seqan-v2.4.0/apps/mason2"
ENTRYPOINT []
