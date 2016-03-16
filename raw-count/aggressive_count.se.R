#!/usr/bin/env Rscript

# Date: 2016-03-17
# calculate raw-count value from bam file.
# If a read matched to multiple genes, then all genes count

args=commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    # files should be in absolute path
    stop("Usage: $0 GTF_FILE BAM_FILE")
}

library(GenomicAlignments)
library(GenomicFeatures)

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

bam <- as(readGAlignments(BamFile(bam.file)), 'GRanges')
res = countOverlaps(exons.by.gene, bam, ignore.strand = TRUE)
write.table(data.frame(names(res), unname(res)),
           row.names=FALSE,
           col.names=FALSE,
           sep='\t')
