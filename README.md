# Ribo-seq Analysis Pipeline — GitHub-ready repository

This repository contains a **step-by-step, reproducible Ribo-seq analysis pipeline** for the GEO datasets **GSE79864** and **GSE66929**. It includes scripts for:

* extracting SRR accessions from GEO using Biopython
* downloading raw SRA fastq files (SRA Toolkit)
* QC with FastQC / MultiQC
* adapter trimming (cutadapt) and rRNA removal (Bowtie2)
* alignment with STAR
* counting footprints with featureCounts
* P-site offset determination with RiboWaltz (R)
* differential translation (DESeq2 + TE calculations)

---

## Repository structure

```
riboseq-pipeline/
├── README.md
├── environment.yml
├── .gitignore
├── LICENSE
├── scripts/
│   ├── get_srr.py
│   ├── download_sra.sh
│   ├── qc_multiqc.sh
│   ├── trim_and_rrna.sh
│   ├── star_index_and_align.sh
│   └── counts_featurecounts.sh
├── R/
│   ├── psite_analysis.R
│   └── deseq_te.R
└── examples/
    └── sample_config.yaml
```

---

> **How to use:** clone the repo, create the conda environment from `environment.yml`, edit `examples/sample_config.yaml` to point to your SRR list / genome files, then run the shell scripts in order.

---

## File contents

---

### `README.md`

````markdown
# Ribo-seq Pipeline (GSE79864 & GSE66929)

## Overview

This repo provides everything needed to go from GEO -> SRR -> FASTQ -> QC -> trimmed/rRNA-depleted reads -> STAR alignment -> counts -> P-site analysis -> differential translation.

## Quickstart

1. Create conda environment:

```bash
conda env create -f environment.yml
conda activate ribo-seq
````

2. Extract SRR list (replace email):

```bash
python3 scripts/get_srr.py --geo GSE79864 GSE66929 --email you@example.com --out srr_list.txt
```

3. Download raw data (SRA Toolkit must be installed):

```bash
bash scripts/download_sra.sh srr_list.txt
```

4. Run QC:

```bash
bash scripts/qc_multiqc.sh ./*.fastq.gz
```

5. Trim adapters and remove rRNA:

```bash
bash scripts/trim_and_rrna.sh srr_list.txt path/to/rRNA.fa
```

6. Build STAR index and align:

```bash
bash scripts/star_index_and_align.sh /path/to/genome.fa /path/to/annotation.gtf srr_sample.fastq.gz
```

7. Generate counts:

```bash
bash scripts/counts_featurecounts.sh /path/to/annotation.gtf aligned.bam
```

8. P-site analysis & plotting (R):

```bash
Rscript R/psite_analysis.R aligned.bam /path/to/annotation.gtf
```

9. Differential translation (R):

```bash
Rscript R/deseq_te.R counts_matrix.csv rna_counts_matrix.csv sample_metadata.csv
```

## Notes

* For large datasets use a cloud VM or HPC.
* Adjust parameters in scripts for single-end vs paired-end data.

```
```

---

### `environment.yml`

```yaml
name: ribo-seq
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - biopython
  - sra-tools
  - fastqc
  - multiqc
  - cutadapt
  - bowtie2
  - star
  - subread  # provides featureCounts
  - r-base
  - r-essentials
  - r-ribowaltz
  - r-deseq2
  - r-pheatmap
  - r-optparse
  - samtools
  - wget
  - gzip
  - pip
  - pip:
    - pyyaml
```

---

### `.gitignore`

```gitignore
.env*
__pycache__/
*.pyc
*.sra
*.fastq
*.fastq.gz
qc_reports/
star_index/
*.bam
*.bai
counts.txt
```

---

### `LICENSE`

```text
MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction... (standard MIT text)
```

---

### `scripts/get_srr.py`

```python
#!/usr/bin/env python3
"""
get_srr.py

Uses Biopython Entrez to retrieve SRR run accessions for given GEO GSE IDs.
Outputs a simple list (one SRR per line).
"""
import argparse
from Bio import Entrez
import time


def gse_to_srr(gse, email):
    Entrez.email = email
    # search GEO GDS db for the GSE
    handle = Entrez.esearch(db="gds", term=gse)
    rec = Entrez.read(handle)
    handle.close()
    if not rec['IdList']:
        return []
    gds_id = rec['IdList'][0]
    # link to sra
    handle = Entrez.elink(dbfrom='gds', db='sra', id=gds_id)
    linkrec = Entrez.read(handle)
    handle.close()
    srrs = []
    if linkrec and linkrec[0].get('LinkSetDb'):
        links = linkrec[0]['LinkSetDb'][0]['Link']
        for link in links:
            sra_id = link['Id']
            # fetch runinfo
            ef = Entrez.efetch(db='sra', id=sra_id, rettype='runinfo', retmode='text')
            runinfo = ef.read()
            ef.close()
            # parse Run column lines
            for line in runinfo.splitlines():
                if line.startswith('Run,'):
                    continue
                parts = line.split(',')
                if parts:
                    run = parts[0]
                    if run.startswith('SRR'):
                        srrs.append(run)
            time.sleep(0.4)
    return sorted(list(set(srrs)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--geo', nargs='+', required=True, help='GEO GSE ids')
    p.add_argument('--email', required=True, help='Entrez email')
    p.add_argument('--out', required=True, help='output file for SRR list')
    args = p.parse_args()

    all_srr = []
    for g in args.geo:
        print(f'Fetching SRR for {g}...')
        s = gse_to_srr(g, args.email)
        print(f'  found {len(s)} runs')
        all_srr.extend(s)

    all_srr = sorted(list(set(all_srr)))
    with open(args.out, 'w') as fh:
        for r in all_srr:
            fh.write(r + '\n')
    print(f'Wrote {len(all_srr)} SRR ids to {args.out}')
```

---

### `scripts/download_sra.sh`

```bash
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
```

> Note: `--split-3` handles single/paired, writes single file for single-end and two for paired (you may prefer `--split-files`).

---

### `scripts/qc_multiqc.sh`

```bash
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
```

---

### `scripts/trim_and_rrna.sh`

```bash
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
```

---

### `scripts/star_index_and_align.sh`

```bash
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
```

---

### `scripts/counts_featurecounts.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

GTF=${1}
shift
BAMS="$@"

featureCounts -T 8 -t CDS -g gene_id -a $GTF -o counts.txt $BAMS

# output: counts.txt
```

---

### `R/psite_analysis.R`

```r
#!/usr/bin/env Rscript
# psite_analysis.R
# Usage: Rscript R/psite_analysis.R sample.bam annotation.gtf

suppressPackageStartupMessages({
  library(riboWaltz)
  library(optparse)
})

option_list = list(
  make_option(c("-b", "--bam"), type="character", help="Aligned BAM file"),
  make_option(c("-g", "--gtf"), type="character", help="GTF annotation file")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

bam = opt$bam
gtf = opt$gtf

# read BAM into riboWaltz
reads_list = bamtolist(bam.files = bam, transcript.annotation = gtf)
ps = psite(reads_list)

# plot per-length P-site distribution
pdf('psite_length_distribution.pdf')
plot_psite(ps)
dev.off()

# metagene plot
pdf('metagene_start_stop.pdf')
plot_metagene(reads_list, lengths = 25:35)
dev.off()

cat('P-site analysis completed. Plots written to current directory.\n')
```

---

### `R/deseq_te.R`

```r
#!/usr/bin/env Rscript
# deseq_te.R
# Usage: Rscript R/deseq_te.R ribo_counts.csv rna_counts.csv sample_metadata.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

option_list = list(
  make_option(c("-ribo"), type="character", help="Ribo counts CSV"),
  make_option(c("-rna"), type="character", help="RNA counts CSV"),
  make_option(c("-meta"), type="character", help="Sample metadata CSV")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ribo_file = opt$ribo
rna_file = opt$rna
meta_file = opt$meta

ribo = read.csv(ribo_file, row.names=1)
rna = read.csv(rna_file, row.names=1)
meta = read.csv(meta_file, row.names=1)

# Align columns by samples (user must ensure samples match)
common_samples = intersect(colnames(ribo), colnames(rna))
ribo = ribo[, common_samples]
rna = rna[, common_samples]

# DESeq2 on Ribo
dds_ribo = DESeqDataSetFromMatrix(countData = ribo, colData = meta[common_samples, ], design = ~ condition)
dds_ribo = DESeq(dds_ribo)
res_ribo = results(dds_ribo)

# DESeq2 on RNA
dds_rna = DESeqDataSetFromMatrix(countData = rna, colData = meta[common_samples, ], design = ~ condition)
dds_rna = DESeq(dds_rna)
res_rna = results(dds_rna)

# Translational efficiency (log2 Ribo / RNA)
library(matrixStats)
log2_ribo = log2(as.matrix(ribo) + 1)
log2_rna = log2(as.matrix(rna) + 1)
te = log2_ribo - log2_rna

# for each gene, do a simple t-test (quick heuristic)
pvals = apply(te, 1, function(x) tryCatch(t.test(x ~ meta[common_samples, 'condition'])$p.value, error=function(e) NA))
adj = p.adjust(pvals, method='BH')
res_te = data.frame(gene=rownames(te), pvalue=pvals, padj=adj)

write.csv(as.data.frame(res_ribo), file='ribo_deseq_results.csv')
write.csv(as.data.frame(res_rna), file='rna_deseq_results.csv')
write.csv(res_te, file='te_results.csv')

cat('Differential translation and TE results written: ribo_deseq_results.csv, rna_deseq_results.csv, te_results.csv\n')
```

---

### `examples/sample_config.yaml`

```yaml
srr_list: srr_list.txt
rRNA_fasta: path/to/rRNA.fa
genome_fasta: path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa
annotation_gtf: path/to/Homo_sapiens.GRCh38.109.gtf
threads: 8
```

---

## Final notes & next steps

* The scripts are intentionally opinionated and conservative (single-ended flow, example adapter). Edit them to match paired-end reads or library-specific adapters.
* Consider replacing `fastq-dump` with `fasterq-dump` and using `--split-files` for paired-end datasets.
* For large-scale runs, run downloads and STAR on cloud VMs or an HPC scheduler (SLURM). The `environment.yml` will help reproducibility.

---

If you'd like, I can:

* produce these files as individual files in a zip for direct download, or
* push them to a GitHub repo for you (if you provide a repo name and a GitHub token), or
* generate a Colab notebook that demonstrates the early steps interactively.
