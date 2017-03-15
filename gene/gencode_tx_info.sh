#!/usr/bin/env bash

# Date: 2017-03-15
# extract the transcript and gene information from Gencode annotation file
# write to stdout

set -ue
gencode=${1:-null}
if [ "$gencode" == "null" ]; then
    echo "Usage: $0 <gencode_annotation_gtf>"
    exit 1
fi

echo -e "gene_id\ttranscript_id\tgene_type\tgene_status\tgene_name\ttranscript_type\ttranscript_status\ttranscript_name"
zcat $gencode \
    | awk -F'\t' '$3=="transcript"{print $9}' \
    | cut -d';' -f 1-8 \
    | grep -oP '(?<=")[^"]*"' \
    | sed 's/"$//' \
    | paste - - - - - - - -
