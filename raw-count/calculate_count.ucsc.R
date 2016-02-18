#!/usr/bin/env Rscript

# Date: 2016-02-18
# calculate raw-count value from bam file
# using TxDb.Hsapiens.UCSC.hg38.knownGene

args=commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    # files should be in absolute path
    stop("Usage: $0 BAM_FILE")
}

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

bam.file <- args[1]
stopifnot(file.exists(bam.file))

exons.by.gene <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")

se <- summarizeOverlaps(exons.by.gene, BamFile(bam.file),
			mode="Union",
			singleEnd=FALSE,
			ignore.strand=TRUE,
			fragments=TRUE)

write.csv(assay(se), col.names=FALSE, quote=FALSE, sep='\t')
