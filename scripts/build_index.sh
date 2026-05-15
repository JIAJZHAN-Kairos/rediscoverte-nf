#!/usr/bin/env bash
# Build REdiscoverTE Salmon index on a beefy box (EC2 c5.4xlarge+) once.
# Uploads to S3 for use as --salmon_index in the Nextflow pipeline.
#
# Prerequisites:
#   - salmon >= 1.10 (conda: mamba install -c bioconda salmon)
#   - aws cli configured
#   - ~120 GB free disk
#   - REdiscoverTE_whole_transcriptome_hg38.fa.xz from Genentech
#     (in rollup_annotation/ of REdiscoverTE_1.0.1.tar.gz)
#
# Usage:
#   ./build_index.sh <path_to_fa.xz> <s3_destination_dir>
# Example:
#   ./build_index.sh \
#       /data/REdiscoverTE_whole_transcriptome_hg38.fa.xz \
#       s3://grimmond-research-nextflow-980504796380-ap-southeast-2-an/refs/REdiscoverTE/salmon_index/

set -euo pipefail

FA_XZ="${1:?usage: $0 <FASTA.xz> <S3_DEST>}"
S3_DEST="${2:?usage: $0 <FASTA.xz> <S3_DEST>}"
THREADS="${THREADS:-$(nproc)}"
WORK="${WORK:-/tmp/rediscoverte_index_build}"

mkdir -p "$WORK"
cd "$WORK"

echo "[1/4] Decompressing FASTA ..."
xz --keep --decompress --threads=0 -c "$FA_XZ" > REdiscoverTE_whole_transcriptome_hg38.fa
ls -lh REdiscoverTE_whole_transcriptome_hg38.fa

echo "[2/4] Building Salmon index (this is the slow part) ..."
# -k 31           : default kmer for reads >75bp
# --keepDuplicates: CRITICAL for TE quantification; keeps near-identical RE copies
# -p $THREADS     : parallel
salmon index \
    -t REdiscoverTE_whole_transcriptome_hg38.fa \
    -i REdiscoverTE_salmon_index \
    -k 31 \
    --keepDuplicates \
    -p "$THREADS"

echo "[3/4] Index built. Size:"
du -sh REdiscoverTE_salmon_index

echo "[4/4] Uploading to S3: $S3_DEST"
aws s3 cp --recursive REdiscoverTE_salmon_index/ "$S3_DEST"

echo
echo "Done. Use this path in nextflow.config / Seqera launch:"
echo "    --salmon_index $S3_DEST"
