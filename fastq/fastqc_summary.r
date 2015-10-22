#!/usr/bin/env Rscript

# Date  : 2015-10-22
# 汇总 fastqc 的输出报告
args=commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
    cat("Usage: $0 [summary.txt]...\n")
    q()
}

library(plyr, quietly=TRUE, warn.conflicts=FALSE)
library(reshape2, quietly=TRUE, warn.conflicts=FALSE)

# 读取所有的 summary.txt 文件
summary <- ldply(args, read.table, sep="\t", stringsAsFactors=FALSE)
summary.attr <- unique(summary[,2])

# 用 cast 方法揉数据
summary <- dcast(summary, V3~V2, value.var="V1")
names(summary)[1] = "file"

# 更改列顺序
summary <- summary[c("file", summary.attr)]
print(summary)

detach('package:plyr')
detach('package:reshape2')
