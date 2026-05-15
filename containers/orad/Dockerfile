# orad container for REdiscoverTE-nf
# Bundles Illumina DRAGEN ORA Decompression v2.7.0 + reference data
# Build:  docker build -t ghcr.io/JIAJZHAN-Kairos/rediscoverte-orad:2.7.0 .
# Push:   docker push    ghcr.io/JIAJZHAN-Kairos/rediscoverte-orad:2.7.0

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ca-certificates gzip tar \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget -q https://webdata.illumina.com/downloads/software/dragen-decompression/orad.2.7.0.linux.tar.gz \
    && tar xzf orad.2.7.0.linux.tar.gz \
    && rm orad.2.7.0.linux.tar.gz \
    && chmod +x /opt/orad.2.7.0.linux/orad

# orad searches $ORA_REF_PATH or $HOME/oradata for refbin
ENV PATH="/opt/orad.2.7.0.linux:${PATH}"
ENV ORA_REF_PATH="/opt/orad.2.7.0.linux/oradata"

# sanity check
RUN orad -D || true

ENTRYPOINT []
CMD ["bash"]
