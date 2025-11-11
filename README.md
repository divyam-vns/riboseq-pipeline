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

---

## Final notes & next steps

* The scripts are intentionally opinionated and conservative (single-ended flow, example adapter). Edit them to match paired-end reads or library-specific adapters.
* Consider replacing `fastq-dump` with `fasterq-dump` and using `--split-files` for paired-end datasets.
* For large-scale runs, run downloads and STAR on cloud VMs or an HPC scheduler (SLURM). The `environment.yml` will help reproducibility.

---
