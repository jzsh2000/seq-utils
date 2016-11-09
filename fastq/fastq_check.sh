#!/bin/bash

# Author: JIN Xiaoyang
# Date  : 2015-10-01
# 检查作为参数的 fq.gz 文件并给出简单统计信息
# 若输出结果中以红字显示 ERROR，说明 gzip 压缩文件出错

if [ $# -eq 0 ] || [ "$1" == "-h" ] ; then
    echo "usage: $0 [fastq_file]..."
    echo "check the compressed file integrity"
    exit 0
fi

for fq_file in $@
do
    [ ! -f $fq_file ]&&continue

    filesize=$(du -sh $fq_file | grep -oP '^[^\s]*')
    if file $fq_file | grep -qiP '\bgzip\b'; then
	filetype="gzip"
	if gzip -t $fq_file 2>/dev/null; then
	    filestat="OK"
	else
	    filestat="[31mERROR[0m"
	fi
    else
	filetype="text"
	filestat="OK"
    fi

    echo -e "$fq_file\t$filetype\t$filesize\t$filestat"
done
