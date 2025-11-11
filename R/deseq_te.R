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
