#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2017-06-24
# Get the dimension info of a (gene expression) matrix file

# TODO: preview first a few lines and fields in file

set -ueo pipefail

file="$1"
file_info=$(file "$file")

# === guess the input file type
file_type='text'
if echo "$file_info" | grep -wFq 'gzip'; then
    file_type='gzip'
elif echo "$file_info" | grep -wFq 'bzip2'; then
    file_type='bzip2'
else
    file_type='text'
fi

# === analyze delimiter
delimiter='comma'
line_1=''
line_2=''
n_row=0
n_col=0
if [ "$file_type" == 'gzip' ]; then
    line_1=$(zcat $file | sed -n '1p')
    line_2=$(zcat $file | sed -n '2p')
    n_row=$(zcat $file | wc -l)
elif [ "$file_type" == 'bzip2' ]; then
    line_1=$(bzcat $file | sed -n '1p')
    line_2=$(bzcat $file | sed -n '2p')
    n_row=$(bzcat $file | wc -l)
else
    line_1=$(cat $file | sed -n '1p')
    line_2=$(cat $file | sed -n '2p')
    n_row=$(cat $file | wc -l)
fi

tab_count_1=$(echo "$line_1" | tr '\t' '\n' | wc -l)
tab_count_2=$(echo "$line_2" | tr '\t' '\n' | wc -l)
com_count_1=$(echo "$line_1" | tr ',' '\n' | wc -l)
com_count_2=$(echo "$line_2" | tr ',' '\n' | wc -l)

if [ $tab_count_1 -eq $tab_count_2 ] && [ $tab_count_1 -gt 1 ]; then
    delimiter='table'
    n_col=$tab_count_1
elif [ $com_count_1 -eq $com_count_2 ] && [ $com_count_1 -gt 1 ]; then
    delimiter='comma'
    n_col=$com_count_1
else
    echo 'Unable to determine field delimiter, abort...'
    exit 1
fi

# === get dimensions
echo "$n_row rows"
echo "$n_col columns"
