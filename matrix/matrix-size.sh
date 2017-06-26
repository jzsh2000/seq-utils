#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2017-06-24
# Get the dimension info of a (gene expression) matrix file

set -ue

file="$1"
file_info=$(file "$file")
file_size=$(ls -lh "$file" | awk '{print $5}')

# === guess the input file type
file_type='text'
if echo "$file_info" | grep -wFq 'gzip'; then
    file_type='gzip'
    # echo 'Detect gzip file' >&2
elif echo "$file_info" | grep -wFq 'bzip2'; then
    file_type='bzip2'
    # echo 'Detect bzip2 file' >&2
else
    file_type='text'
    # echo 'Detect plain text file' >&2
fi

# === analyze delimiter
delimiter='comma'
line_1=''
line_2=''
line_top10=''
n_row=0
n_col=0

if [ "$file_type" == 'gzip' ]; then
    line_1=$(zcat $file | sed -n '1p')
    line_2=$(zcat $file | sed -n '2p')
    line_top10=$(zcat $file | head)
    n_row=$(zcat $file | wc -l)
elif [ "$file_type" == 'bzip2' ]; then
    line_1=$(bzcat $file | sed -n '1p')
    line_2=$(bzcat $file | sed -n '2p')
    line_top10=$(bzcat $file | head)
    n_row=$(bzcat $file | wc -l)
else
    line_1=$(cat $file | sed -n '1p')
    line_2=$(cat $file | sed -n '2p')
    line_top10=$(cat $file | head)
    n_row=$(cat $file | wc -l)
fi

tab_count_1=$(echo "$line_1" | tr '\t' '\n' | wc -l)
tab_count_2=$(echo "$line_2" | tr '\t' '\n' | wc -l)
com_count_1=$(echo "$line_1" | tr ',' '\n' | wc -l)
com_count_2=$(echo "$line_2" | tr ',' '\n' | wc -l)

snippet=''
if [ $tab_count_1 -eq $tab_count_2 ] && [ $tab_count_1 -gt 1 ]; then
    delimiter='table'
    snippet=$(echo "$line_top10" | cut -f1-5 | csvlook -t)
    n_col=$tab_count_1
elif [ $com_count_1 -eq $com_count_2 ] && [ $com_count_1 -gt 1 ]; then
    delimiter='comma'
    snippet=$(echo "$line_top10" | cut -d',' -f1-5 | csvlook)
    n_col=$com_count_1
else
    echo 'Unable to determine field delimiter, abort...'
    exit 1
fi

# === get dimensions
acolor='\033[0;36m'
ncolor='\033[0m'
file_size_color=${acolor}${file_size}${ncolor}
n_row_color=${acolor}${n_row}${ncolor}
n_col_color=${acolor}${n_col}${ncolor}
echo -e "$(basename $file): \
    [ $file_size_color | $n_row_color rows | $n_col_color columns ]"
echo "$snippet"
