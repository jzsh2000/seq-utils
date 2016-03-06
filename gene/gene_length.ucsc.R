#!/usr/bin/env Rscript

# Date: 2016-01-22
# Reference: https://www.biostars.org/p/83901/

# calculate length of each gene
# using TxDb.Hsapiens.UCSC.hg38.knownGene

args=commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    cat("Usage: $0 OUTPUT\n")
    q()
} else {
    output.file = args[1]
}

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene)
exons.list.per.gene <- exonsBy(txdb,by="gene")

gene.symbol = mapIds(org.Hs.eg.db,
    keys=names(exons.list.per.gene),
    keytype="ENTREZID",
    column="SYMBOL")

# gene.exons_num = elementLengths(exons.list.per.gene)
gene.length = sapply(exons.list.per.gene, function(x){sum(width(reduce(x)))})

output <- data.frame(id=names(exons.list.per.gene), symbol=gene.symbol, length=gene.length)
write.csv(output, output.file, row.names=FALSE)
