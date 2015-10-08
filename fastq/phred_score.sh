#! /bin/bash

# 2015-09-29
# guess phred-score of fastq files

if [ $# -eq 0 ]; then
    echo "Usage: $0 [FASTQ_FILE]..."
fi

for file in $@
do
    [ ! -f $file ]&&continue

    echo -n "$file	"
    if gzip -t $file 2>/dev/null; then
	phred33_line=$(zcat $file | awk 'NR%4==0{print}' | grep -c '[!-?]')
    else
	phred33_line=$(awk 'NR%4==0{print}' $file | grep -c '[!-?]')
    fi

    echo -n "$phred33_line	"
    if [ $phred33_line -eq 0 ]; then
	echo "phred-64"
    else
	echo "phred-33"
    fi
done
