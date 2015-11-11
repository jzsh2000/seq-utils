#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2015-11-10
# The encapsulation of `getSubSeq` program

realdir=$(dirname $(realpath $0))
realbin=$realdir/getSubSeq

# ======================= user config =========================== #

## some default genome fasta file paths
hg19_fa='/Users/jinxy/Database/hg19/genome.fa'
hg38_fa='/Users/jinxy/Database/hg38/genome.fa'

# ======================= end user config =========================== #

function usage() {
    echo "Usage: $0 [-r] chr pos1 pos2 [fa_id | fa_file]"
}

function check_num() {
    if ! echo "$1" | egrep -xq '[0-9]+'; then
	echo "Error: not a number - '$1'" >&2
	exit 1
    fi
}

if [ "$1" == "-r" ]; then
    rev=1
    shift
else
    rev=0
fi

if [ $# -lt 3 ]; then
    usage >&2
    exit 1
fi

chr=$1
pos1=$2
pos2=$3
fa_file=$4

check_num $pos1
check_num $pos2

# pos1 shouldn't be greater than pos2, but if pos1 is really greater than pos2,
# then the output would be reverse complementary sequences. The effect is same
# as add `-r` parament. Note that if `-r` parament is added while pos1 is
# greater than pos2, these two decorations would be neutralized.
if [ $pos1 -gt $pos2 ]; then
    rev=$[1-$rev]
    pos1=$3
    pos2=$2
fi

if [ -z $fa_file ]; then
    # use hg38 if `fa_file` not set
    fa_file=$hg38_fa
elif [ "$fa_file" == "hg19" ]; then
    fa_file=$hg19_fa
elif [ "$fa_file" == "hg38" ]; then
    fa_file=$hg38_fa
fi

if [ ! -f $fa_file ] || [ ! -r $fa_file ]; then
    echo "Error: cannot read file '$fa_file'" >&2
fi

if [ $rev -eq 0 ]; then
    $realbin -q ${chr}:${pos1}-${pos2} -f $fa_file | tr 'a-z' 'A-Z' | fold -w 50
else
    $realbin -q ${chr}:${pos1}-${pos2} -f $fa_file | rev | tr 'a-z' 'A-Z' | tr 'AGCT' 'TCGA' | fold -w 50
fi
