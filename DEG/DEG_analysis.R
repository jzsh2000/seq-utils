#!/usr/bin/env Rscript

# Get differential expression genes according to fpkm or raw-count

# Two pipelines are used:
# tophat => cufflinks => cuffdiff
# tophat => htseq-count => edgeR

# Author: JIN Xiaoyang
# Date  : 2015-10

# get arguments
args=commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
    cat("Usage: $0 [cuffnorm_out_dir] [cuffdiff_out_dir]\n")
    cat("    or $0 [htseq-count_out_dir]\n")
    q()
}

for (arg in args) {
    # get cuffnorm directory
    if(file.exists(file.path(arg, "genes.attr_table")) && file.exists(file.path(arg, "genes.fpkm_table"))) {
	cuffnorm.dir=arg
    }
    # get cuffdiff directory
    if(file.exists(file.path(arg, "gene_exp.diff"))) {
	cuffdiff.dir=arg
    }
    # get HTSeq-count directory
    if(file.exists(file.path(arg, "raw-count.txt"))) {
	htseq.dir=arg
    }
}

if (exists("cuffnorm.dir") && exists("cuffdiff.dir")) {

    library(plyr, quietly=TRUE, warn.conflicts=FALSE)
    library(ggplot2, quietly=TRUE, warn.conflicts=FALSE)

    cuffnorm.attr <- read.delim(file.path(cuffnorm.dir, "genes.attr_table"), stringsAsFactors=FALSE, check.names=FALSE)
    cuffnorm.fpkm <- read.delim(file.path(cuffnorm.dir, "genes.fpkm_table"), stringsAsFactors=FALSE, check.names=FALSE)
    cuffdiff.fpkm <- read.delim(file.path(cuffdiff.dir, "gene_exp.diff"), stringsAsFactors=FALSE, check.names=FALSE)
    cuffdiff.data <- merge(cuffdiff.fpkm, cuffnorm.fpkm, by.x="gene_id", by.y="tracking_id")

    # cuffdiff.data <- subset(cuffdiff.data, gene!="-"&apply(cuffdiff.data[, c("value_1", "value_2")], 1, max)>1)
    # TODO: update filter rule
    cuffdiff.data <- subset(cuffdiff.data, gene!="-"&apply(cuffdiff.data[, 15:ncol(cuffdiff.data)], 1, function(x){sum(x>=1)})>=3)

    # gene expression change - none/upgrade/downgrade
    cuffdiff.data$change = "none"
    cuffdiff.data$change[cuffdiff.data$significant=="yes"&cuffdiff.data$`log2(fold_change)`>0]="upgrade"
    cuffdiff.data$change[cuffdiff.data$significant=="yes"&cuffdiff.data$`log2(fold_change)`<0]="downgrade"

    # cuffdiff DEG list
    cuffdiff.deg <- subset(cuffdiff.data, significant=="yes")

    # calculate DEG numbers per sample pair
    cuffdiff.count = ddply(cuffdiff.deg, .(sample_1, sample_2), nrow)
    names(cuffdiff.count)[3]="deg.count"

    # transform data.frame to matrix
    sample.name=unique(c(cuffdiff.deg$sample_1, cuffdiff.deg$sample_2))
    cuffdiff.count.matrix <- matrix(0, length(sample.name), length(sample.name))
    rownames(cuffdiff.count.matrix)=sample.name
    colnames(cuffdiff.count.matrix)=sample.name

    for (i in 1:nrow(cuffdiff.count)) {
	cuffdiff.count.matrix[cuffdiff.count[i,1], cuffdiff.count[i,2]]=cuffdiff.count[i,3]
	cuffdiff.count.matrix[cuffdiff.count[i,2], cuffdiff.count[i,1]]=cuffdiff.count[i,3]
    }
    cuffdiff.count.matrix=cuffdiff.count.matrix[order(cuffdiff.count.matrix[,1]) ,]
    cuffdiff.count.matrix=cuffdiff.count.matrix[, order(cuffdiff.count.matrix[1,])]

    cat("==[31mcuffdiff_DEG_count[0m==\n")
    print(cuffdiff.count.matrix)

    if(file.exists(file.path(cuffdiff.dir, "gene.label"))) {
	gene.label <- readLines(file.path(cuffdiff.dir, "gene.label"))
	cuffdiff.data$label=""
	cuffdiff.data$label.bool=FALSE
	for (mylabel in gene.label) {
	    # print(which(grepl(paste("\\b\\Q", label, "\\E\\b", sep=""), cuffdiff.data$gene, ignore.case=TRUE, perl=TRUE)))
	    mylabel.line=grepl(paste("\\b\\Q", mylabel, "\\E\\b", sep=""), cuffdiff.data$gene, ignore.case=TRUE, perl=TRUE)
	    cuffdiff.data$label[mylabel.line]=mylabel
	    cuffdiff.data$label.bool[mylabel.line]=TRUE
	}

	# fpkm.plotobj
	fpkm.plotobj=ggplot(cuffdiff.data, aes(x=log2(value_1), y=log2(value_2), color=change, shape=label.bool))+ scale_color_manual(values=c("blue", "gray", "red"))+ scale_shape_manual(values=c(20, 8), guide=FALSE)+ geom_point(alpha=0.8)+ xlab(paste("log2(", cuffdiff.data[, c("sample_1")] ,".fpkm)"))+ ylab(paste("log2(", cuffdiff.data[, c("sample_2")] ,".fpkm)"))

	fpkm.plotobj=fpkm.plotobj+geom_text(aes(x=log2(value_1), y=log2(value_2), label=cuffdiff.data$label), size=1.5, hjust=0, vjust=1, colour="black")+ guides(label.bool=FALSE)
    }
    else {
	# fpkm.plotobj
	fpkm.plotobj=ggplot(cuffdiff.data, aes(x=log2(value_1), y=log2(value_2), color=change))+ scale_color_manual(values=c("blue", "gray", "red"))+ geom_point(alpha=0.8)+ xlab(paste("log2(", cuffdiff.data[, c("sample_1")] ,".fpkm)"))+ ylab(paste("log2(", cuffdiff.data[, c("sample_2")] ,".fpkm)"))
    }

    ggsave(plot = fpkm.plotobj, filename="DEG.fpkm.pdf", height=6, width=8)

}

if (exists("htseq.dir")) {
    readcount.raw <- read.delim(file.path(htseq.dir, "raw-count.txt"), row.names=1, stringsAsFactors=FALSE)
    readcount.raw <- readcount.raw[1:(nrow(readcount.raw)-5),]

    library(ggplot2, quietly=TRUE, warn.conflicts=FALSE)
    library(edgeR, quietly=TRUE, warn.conflicts=FALSE)

    # tail(readcount.raw)
    readcount.sum=colSums(readcount.raw)
    readcount.data <- data.frame(name=names(readcount.sum), reads=unname(readcount.sum), sample=str_extract(names(readcount.sum), regex(".*(?=_)")), rep=str_extract(names(readcount.sum), regex("[^_]*$")))
    print(readcount.data)
    quick.barplot <- ggplot(readcount.data, aes(x=rep, y=reads, fill=sample))+geom_bar(stat="identity", position="dodge")
    ggsave(plot = quick.barplot, filename="reads.count.pdf", height=6, width=8)

    # run edgeR
    edgeR.group=factor(readcount.data$sample)
    sample.name=levels(edgeR.group)

    edgeR.count.matrix <- matrix(0, length(sample.name), length(sample.name))
    rownames(edgeR.count.matrix)=sample.name
    colnames(edgeR.count.matrix)=sample.name

    # edgeR - classic edgeR analysis
    edgeR.y <- DGEList(counts = readcount.raw, group = edgeR.group)
    edgeR.y <- calcNormFactors(edgeR.y)
    edgeR.y <- estimateCommonDisp(edgeR.y)
    edgeR.y <- estimateTagwiseDisp(edgeR.y)

    for (sam1 in sample.name) {
	for (sam2 in sample.name) {
	    if(sam1==sam2) {
		next
	    }
	    edgeR.et <- exactTest(edgeR.y, pair = c(sam1, sam2))
	    edgeR.tt <- topTags(edgeR.et, n=nrow(readcount.raw))
	    edgeR.deg <- rownames(edgeR.tt)[edgeR.tt$table$FDR < 0.05]
	    edgeR.count.matrix[sam1, sam2]=length(edgeR.deg)
	}
    }

    edgeR.count.matrix=edgeR.count.matrix[order(edgeR.count.matrix[,1]) ,]
    edgeR.count.matrix=edgeR.count.matrix[, order(edgeR.count.matrix[1,])]

    cat("==[31medgeR_DEG_count.classic[0m==\n")
    print(edgeR.count.matrix)

    # edgeR - glm edgeR analysis
    edgeR.count.matrix <- matrix(0, length(sample.name), length(sample.name))
    rownames(edgeR.count.matrix)=sample.name
    colnames(edgeR.count.matrix)=sample.name

    edgeR.design <- model.matrix(~0+edgeR.group)
    edgeR.y <- DGEList(counts = readcount.raw, group = edgeR.group)
    edgeR.y <- estimateGLMCommonDisp(edgeR.y, edgeR.design)
    edgeR.y <- estimateGLMTrendedDisp(edgeR.y, edgeR.design)
    edgeR.y <- estimateGLMTagwiseDisp(edgeR.y, edgeR.design)
    edgeR.fit <- glmFit(edgeR.y, edgeR.design)

    for (sam1 in sample.name) {
	for (sam2 in sample.name) {
	    if(sam1==sam2) {
		next
	    }
	    # edgeR.contrast <- rep(0, length(sample.name))
	    # edgeR.contrast[match(sam1, sample.name)]=-1
	    # edgeR.contrast[match(sam2, sample.name)]=1
	    # print(edgeR.contrast)
	    edgeR.contrast=makeContrasts(paste("edgeR.group",sam2,"-","edgeR.group",sam1, sep=""), levels=edgeR.design)

	    edgeR.lrt <- glmLRT(edgeR.fit, contrast=edgeR.contrast)
	    edgeR.tt <- topTags(edgeR.lrt, n=nrow(readcount.raw))
	    edgeR.deg <- rownames(edgeR.tt)[edgeR.tt$table$FDR < 0.05]
	    edgeR.count.matrix[sam1, sam2]=length(edgeR.deg)
	}
    }

    edgeR.count.matrix=edgeR.count.matrix[order(edgeR.count.matrix[,1]) ,]
    edgeR.count.matrix=edgeR.count.matrix[, order(edgeR.count.matrix[1,])]

    cat("==[31medgeR_DEG_count.glm[0m==\n")
    print(edgeR.count.matrix)


}
