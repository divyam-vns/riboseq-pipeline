#!/usr/bin/env bash
set -euo pipefail

GENOME_FA=$1
GTF=$2
THREADS=${3:-8}
GENOME_DIR=${4:-star_index}

mkdir -p $GENOME_DIR

# Build index (run once)
STAR --runThreadN $THREADS --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $GENOME_FA \
     --sjdbGTFfile $GTF \
     --sjdbOverhang 29

# Align all *_no_rRNA.fastq.gz files in cwd
for fq in *_no_rRNA.fastq.gz; do
  sample=$(basename $fq _no_rRNA.fastq.gz)
  echo "Aligning $sample"
  STAR --runThreadN $THREADS \
       --genomeDir $GENOME_DIR \
       --readFilesIn $fq \
       --readFilesCommand zcat \
       --outFileNamePrefix ${sample}_ \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMunmapped Within
  mv ${sample}_Aligned.sortedByCoord.out.bam ${sample}.bam
done
