#!/bin/bash

# Author: JIN Xiaoyang
# Date  : 2015-10-09
# 计算文件中的各字符数，输出到 STDOUT

awk '
{
    for (i=1; i<=length; i++)
	arr[substr($0, i, 1)]++
    }
    END {
    for (i in arr) {
	print i, arr[i]
    }
}' $@
