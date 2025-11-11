#!/usr/bin/env bash
set -euo pipefail

SRR_LIST=${1:-srr_list.txt}
RRNA_FASTA=${2:-rRNA.fa}
THREADS=${3:-4}

# build bowtie2 index for rRNA if not present
if [ ! -f rRNA_index.1.bt2 ]; then
  bowtie2-build $RRNA_FASTA rRNA_index
fi

while read -r srr; do
  echo "Processing $srr"
  IN=${srr}.fastq.gz
  TRIM=${srr}_trim.fastq.gz
  NO_RRNA=${srr}_no_rRNA.fastq.gz

  # adapter trimming (example adapter - adjust as needed)
  cutadapt -a AGATCGGAAGAGC -q 20 -m 25 -o $TRIM $IN

  # filter rRNA: keep reads that DO NOT map to rRNA
  bowtie2 -x rRNA_index -U $TRIM --very-sensitive-local --un $NO_RRNA -S /dev/null -p $THREADS

  echo "Produced: $NO_RRNA"

done < $SRR_LIST
