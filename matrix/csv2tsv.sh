#!/usr/bin/env bash

# Author: JIN Xiaoyang
# Date  : 2017-06-28
# Convert csv file to tsv file

set -ueo pipefail

cat $1 \
    | tr ',' '\t' \
    > ${1%.csv}.tsv