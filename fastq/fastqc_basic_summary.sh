#!/bin/bash

# Date  : 2015-10-22
# 汇总 fastqc 输出的·Basic Statistics·信息 (文件名为 fastqc_data.txt)

## fastqc_data.txt 文件头示例:

## ##FastQC	0.11.3
## >>Basic Statistics	pass
## #Measure	Value
## Filename	DMSO_R1_P1.fq.gz
## File type	Conventional base calls
## Encoding	Illumina 1.5
## Total Sequences	24646662
## Sequences flagged as poor quality	0
## Sequence length	90
## %GC	48

output=""

for fastqc_data in $@
do
    [ ! -f "$fastqc_data" ]&&continue

    # 为了解决换行的问题，不同行先用'@'符号连接，最后再统一换成'\n'
    if [ -z "$output" ]; then
	output="$(grep -m1 '>>' -A8 $fastqc_data | cut -f1 | tail -n+3 | tr '\n' '\t')"
	output="${output}@"
    fi

    output="${output}$(grep -m1 '>>' -A8 $fastqc_data | cut -f2 | tail -n+3 | tr '\n' '\t')@"
done

echo "$output" | tr '@' '\n' | column -s $'\t' -t
