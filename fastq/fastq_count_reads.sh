#!/bin/bash

# Author: JIN Xiaoyang
# Date  : 2015-10-12
# counts the number of reads in a fastq file

if [ $# -eq 0 ] ; then
    echo "usage: $0 [fastq_file]..."
    echo "counts the number of reads in a fastq file"
    exit 0
fi
 
for fq_file in $@
do
    if [ -z "$fq_file" ]||[ ! -f "$fq_file" ]; then
	continue
    fi

    if file "$fq_file" | grep -qiP '\bgzip\b'; then
	lines=$(zcat $fq_file | wc -l | cut -d' ' -f1)
    else
	lines=$(wc -l $fq_file | cut -d' ' -f1)
    fi

    count=$[$lines / 4]
    
    if [ $# -eq 1 ]; then
	echo $count
    else
	printf "%-20d%s\n" $count $fq_file
    fi
done
