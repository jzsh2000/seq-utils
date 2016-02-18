#!/usr/bin/env Rscript

# Date: 2016-01-22
# Reference: https://www.biostars.org/p/83901/

# calculate length of each gene, assume gtf or gff3 file as input
# for human annotation, this script may run a long time
args=commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("Usage: $0 GTF_FILE OUTPUT\n")
    q()
}

library(GenomicFeatures, quietly=TRUE, warn.conflicts=FALSE)
txdb <- makeTxDbFromGFF(args[1])
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

write.table(as.matrix(exonic.gene.sizes), args[2], sep="\t", quote=F, col.names=F)
