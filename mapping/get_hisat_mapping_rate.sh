#!/usr/bin/env bash

# Date: 2016-01-17
# Read mapping rate from hisat output file

for file in $@
do
    [ ! -f $file ] && continue

    echo -ne "$(basename $file)\t"
    tail -1 $file | grep -o '^[^ ]*'
done
