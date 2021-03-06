#!/bin/bash

# Author: JIN Xiaoyang
# Date  : 2015-09-29
# merge the output files of htseq-count (with suffix *.count)

file_list=($@)
name_list=(${file_list[@]/*\//})
name_list=(${name_list[@]/.count/})

echo -e "\t${name_list[@]}" | tr ' ' '\t'
paste ${file_list[@]} | awk '{printf("%s",$1); for(i=2;i<=NF;i=i+2){ printf("\t%s",$i)} printf("%c","\n")}'
