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
