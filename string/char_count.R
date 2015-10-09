#!/usr/bin/env Rscript

# Author: JIN Xiaoyang
# Date  : 2015-10-09
# 计算文件中的各字符数，输出到 STDOUT
# 输出结果相对紧凑，不能直接读取 STDIN

# get arguments
args=commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
    cat("Usage: $0 FILE...\n")
    q()
}

library(stringr, quietly=TRUE, warn.conflicts=FALSE)

for (arg in args) {
    if(file.exists(arg)) {
	file.content <- readLines(arg)

	if(length(args)>1) {
	    cat(paste("[1m>", arg, "[0m"))
	}
	print(table(unlist(str_extract_all(file.content,regex(".")))))
    }
}

detach('package:stringr')
