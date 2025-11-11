#!/usr/bin/env bash
set -euo pipefail

SRR_LIST=${1:-srr_list.txt}
THREADS=${2:-4}

while read -r srr; do
  echo "Downloading $srr"
  # use prefetch to get .sra file
  prefetch --max-size 100G $srr
  # convert to fastq.gz
  fastq-dump --split-3 --gzip --stdout $srr > ${srr}.fastq.gz || true
  # optionally move sra to archive
  if [ -f ${srr}.sra ]; then
    mv ${srr}.sra sra_arch/
  fi
done < $SRR_LIST

echo "All downloads complete."

## Note: --split-3 handles single/paired, writes single file for single-end and two for paired (you may prefer --split-files).
