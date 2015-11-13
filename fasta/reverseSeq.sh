#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2015-11-12
#
# reverse genome fasta file
#=========================#
# Depends: GNU coreutils  #
#   * csplit              #
#   * head                #
#   * tail                #
#   * rev                 #
#   * tac                 #
#   * tr                  #
#   * fold                #
#   * rm                  #
#=========================#

function usage() {
    echo "Usage: $0 <fa_file> [width]"
}

if [ $# -eq 0 ]; then
    usage
    exit 0
fi

fa_file=$1
if [ ! -f $fa_file ]; then
    echo "Error: cannot open file '$fa_file'" >&2
    exit 1
fi

width=$2
if [ -z "$width" ]; then
    width=50
fi

csplit -b "%05d" -f "/tmp/fa.tmp.$$." -z -s $fa_file '/^>/' '{*}'

for fa_tmp_file in $(ls /tmp/fa.tmp.$$.?????)
do
    head -1 $fa_tmp_file
    tail -n+2 $fa_tmp_file | rev | tac | tr 'agctAGCT' 'tcgaTCGA' | tr -d '\n' | fold -w $width
    echo
done

rm -f /tmp/fa.tmp.$$.?????
