#!/usr/bin/env Rscript

# Date: 2016-01-28
# calculate raw-count value from bam file

args=commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    # files should be in absolute path
    stop("Usage: $0 GTF_FILE BAM_FILE")
}

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)

gtf.file  <- args[1]
bam.file <- args[2]
tmpfile <- '/tmp/.gtf.rda'

stopifnot(file.exists(bam.file), file.exists(gtf.file))

if(file.exists(tmpfile)) {
    load(tmpfile)
} else {
    txdb <- makeTxDbFromGFF(gtf.file)
    exons.by.gene <- exonsBy(txdb, by="gene")
    save(txdb, exons.by.gene, file = tmpfile)
}

se <- summarizeOverlaps(exons.by.gene, BamFile(bam.file),
			mode="Union",
			singleEnd=FALSE,
			ignore.strand=TRUE,
			fragments=TRUE)

write.table(assay(se), col.names=FALSE, quote=FALSE, sep='\t')
