#!/usr/bin/env Rscript

# Date: 2016-03-16
# calculate raw-count value from bam file

args=commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    # files should be in absolute path
    stop("Usage: $0 GTF_FILE BAM_FILE...")
}

library(Rsubread)

gtf.file  <- args[1]
bam.files <- args[-1]

stopifnot(file.exists(bam.files), file.exists(gtf.file))

fc <- featureCounts(bam.files, annot.ext = gtf.file,
    isGTFAnnotationFile = TRUE, 
    countMultiMappingReads = TRUE)

write.table(fc$counts, col.names=FALSE, sep='\t')
