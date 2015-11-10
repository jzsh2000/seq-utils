#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2015-11-10
# The encapsulation of `getSubSeq` program

realdir=$(dirname $(realpath $0))
realbin=$realdir/getSubSeq

# ======================= user config =========================== #

## some default genome fasta file paths
hg19_fa='/Volumes/STDT4000300/Database/iGenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
hg38_fa='/Volumes/STDT4000300/Database/iGenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'

# ======================= end user config =========================== #

function usage() {
    echo "Usage: $0 chr pos1 pos2 [fa_id | fa_file]"
}

chr=$1
pos1=$2
pos2=$3
fa_file=$4

if [ $# -lt 3 ]; then
    usage
fi

if [ -z $fa_file ] || [ "$fa_file" == "hg19" ]; then
    # use hg19 if `fa_file` not set
    fa_file=$hg19_fa
elif [ "$fa_file" == "hg38" ]; then
    fa_file=$hg38_fa
fi

if [ ! -f $fa_file ] || [ ! -r $fa_file ]; then
    echo "Error: cannot read file '$fa_file'" >&2
fi

$realbin -q ${chr}:${pos1}-${pos2} -f $fa_file | fold
