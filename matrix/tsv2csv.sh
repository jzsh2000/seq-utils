#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2017-06-28
# Convert tsv file to csv file

set -ueo pipefail

cat $1 \
    | tr '\t' ',' \
    > ${1%.tsv}.csv