# Nirvana 3.18.1 (Illumina) on the .NET 6 runtime: VCF annotation for BIAS-2015.
# Annotation data (Cache, SupplementaryAnnotation, References) is mounted, never baked in.
# The .NET 6 image is Debian 11, whose security archive is gone from the mirrors, so the release
# is fetched in a Debian 12 build stage and only copied into the runtime image.
#   docker build -t clinval-nirvana:3.18.1 containers/nirvana
FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587 AS fetch
ARG NIRVANA_VERSION=3.18.1
RUN apt-get update && apt-get install -y --no-install-recommends curl unzip ca-certificates \
    && curl -fsSL -o /tmp/nirvana.zip \
       https://github.com/Illumina/Nirvana/releases/download/v${NIRVANA_VERSION}/Nirvana-${NIRVANA_VERSION}-net6.0.zip \
    && mkdir -p /opt/nirvana && unzip -q /tmp/nirvana.zip -d /opt/nirvana \
    && NIRVANA_DLL=$(find /opt/nirvana -name Nirvana.dll | head -1) \
    && DOWNLOADER_DLL=$(find /opt/nirvana -name Downloader.dll | head -1) \
    && printf '#!/bin/sh\nif [ "$1" = "--version" ]; then echo "Nirvana %s"; exit 0; fi\nexec dotnet %s "$@"\n' \
       "${NIRVANA_VERSION}" "${NIRVANA_DLL}" > /opt/nirvana/nirvana \
    && printf '#!/bin/sh\nexec dotnet %s "$@"\n' "${DOWNLOADER_DLL}" > /opt/nirvana/nirvana_downloader \
    && chmod +x /opt/nirvana/nirvana /opt/nirvana/nirvana_downloader

FROM mcr.microsoft.com/dotnet/runtime:6.0@sha256:75752d6e7d17246d939d8f0717a2daf14ca276e1fd6921caa489a70e24498881
COPY --from=fetch /opt/nirvana /opt/nirvana
RUN ln -s /opt/nirvana/nirvana /usr/local/bin/nirvana \
    && ln -s /opt/nirvana/nirvana_downloader /usr/local/bin/nirvana_downloader
ENV DOTNET_CLI_TELEMETRY_OPTOUT=1 HOME=/tmp
LABEL org.opencontainers.image.title="nirvana" org.opencontainers.image.version="3.18.1" \
      org.opencontainers.image.source="https://github.com/Illumina/Nirvana/releases/tag/v3.18.1"
ENTRYPOINT []
