#!/bin/bash

# 2015-09-25
# 根据 fastqc 输出的 html 文件判断 fastq 文件版本（phred-score）

if [ $# -eq 0 ]; then
    echo "Usage: $0 [HTML]..."
fi

for html in $@ 
do
    [ ! -f "$html" ]&&continue
    echo -ne "${html##*/}\t"
    grep -m1 -oP 'Illumina [.0-9]+' $html
done
