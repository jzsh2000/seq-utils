#!/usr/bin/env Rscript

# Date: 2016-01-22
# Reference: https://www.biostars.org/p/83901/

# calculate length of each gene
# using TxDb.Hsapiens.UCSC.hg38.knownGene

args=commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    warning("write table to STDOUT")
    output.file = ""
} else {
    output.file = args[1]
}

suppressMessages(library(GenomicFeatures))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(stringr))

txdb <- keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb.exons <- exonsBy(txdb,by="gene")

getMapInfo <- function(gene_ids, column) {
    dat = select(org.Hs.eg.db,
                   keys = gene_ids,
                   keytype = "ENTREZID",
                   columns = column)
    tapply(dat[[column]], dat$ENTREZID, function(vec) {
        if(is.na(vec)) {
            NA
        } else {
            paste(vec, collapse = ' | ')
        }
    })
}

geneinfo <- data.frame(
    id = names(txdb.exons),
    symbol = getMapInfo(names(txdb.exons), "SYMBOL"),
    map = getMapInfo(names(txdb.exons), "MAP"),
    alias = getMapInfo(names(txdb.exons), "ALIAS"),
    name = getMapInfo(names(txdb.exons), "GENENAME"),
    # exons = elementLengths(txdb.exons),
    # txs = elementLengths(transcriptsBy(txdb, by="gene")[names(txdb.exons)]),
    length = sapply(txdb.exons, function(x) {sum(width(reduce(x)))})
) 

geneinfo = geneinfo[complete.cases(geneinfo), ]
geneinfo$chr = sapply(str_split(geneinfo$map, ' \\| '), function(strlist) {
    paste(paste0("chr", str_extract(strlist, pattern = "[0-9]+|[XY]")),
        collapse = ' | ')
})

write.csv(geneinfo, output.file, row.names=FALSE)
