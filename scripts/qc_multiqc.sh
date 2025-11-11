#!/usr/bin/env bash
set -euo pipefail

OUTDIR=${1:-qc_reports}
mkdir -p $OUTDIR

for fq in "$@"; do
  if [[ $fq == *.fastq.gz ]]; then
    fastqc -o $OUTDIR $fq || true
  fi
done

multiqc $OUTDIR -o ${OUTDIR}/multiqc_report
