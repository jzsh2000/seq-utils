#!/usr/bin/env Rscript

# Usage: pca.R <count_data> <sample_info> <output_file>
#              [-l <label>] [-c <color>] [-s <shape>]
# Example: ./pca.R test/count_data.csv test/sample_info.csv test.pdf \
#                  -l label
# ========== read parameters
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: pca.R <count_data> <sample_info> <output_file>\n")
    cat("             [-l <label>] [-c <color>] [-s <shape>]\n")
    q(save = "no")
}

count_data <- read.csv(args[1], row.names = 1, stringsAsFactors = FALSE)
sample_info <- read.csv(args[2], stringsAsFactors = FALSE)
output_file = args[3]

args = args[-c(1:3)]
label = NULL
color = NULL
shape = NULL
while (length(args) > 0) {
    if (args[1] == "-l") {
        label = args[2]
        args = args[-c(1:2)]
    } else if (args[1] == "-c") {
        color = args[2]
        args = args[-c(1:2)]
    } else if (args[1] == "-s") {
        shape = args[2]
        sample_info[,shape] = as.character(sample_info[,shape])
        args = args[-c(1:2)]
    } else {
        break
    }
}

cat(paste("color:", color, "\n"))
cat(paste("label:", label, "\n"))
cat(paste("shape:", shape, "\n"))

# ========== START

library(tidyverse)
library(ggrepel)
library(genefilter)
library(DESeq2)

# count_data: raw count matrix, colnames are the sample labels
# sample_info: data frame of sample information
# label: column name in sample_info which represents label in plot
# color: column name in sample_info which represents color in plot
# shape: column name in sample_info which represents shape in plot
count_data = floor(count_data)
sample_info$sample_idx = factor(seq(nrow(sample_info)))
design_formula = "~ sample_idx"

dds <- DESeqDataSetFromMatrix(count_data, colData = sample_info,
                              design = as.formula(design_formula))
rld = rlog(dds)

sample_group = c(label, color, shape)
if (is.null(sample_group) || length(sample_group) == 0) {
    sample_group = "sample_idx"
}
data <- plotPCA(rld,
                intgroup = sample_group,
                returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
p = ggplot(data,
           aes_string(x = "PC1", y = "PC2")) +
    geom_point(aes_string(color = color,
                          shape = shape),
               alpha = .7) +
    ggtitle("PCA") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_bw()

if (!is.null(label)) {
    p = p +
        geom_label_repel(aes_string(x = "PC1", y = "PC2",
                                    label = label),
                         size = 2,
                         point.padding = unit(0.3, "lines"))
}
ggsave(p, filename = output_file)

# ========== END
