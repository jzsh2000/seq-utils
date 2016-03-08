#!/usr/bin/env bash

# Date: 2016-03-08
# Read number of mapped reads from hisat output
# write to stdout as a csv file

echo -e "filename,total_reads,mapped_reads,mapping_rate"

for file in $@
do
    [ ! -f $file ] && continue

    echo -ne "$(basename $file),"

    total_frags=$(head -1 $file | grep -oP '^\d+')
    bad_frags=$(grep 'aligned 0 times$' $file | grep -oPm1 '^ +\d+')

    echo -ne "${total_frags},"

    if [ $(cat $file | wc -l) -gt 10 ]; then
        mode=paired
        echo -ne "$[${total_frags} - ${bad_frags} / 2],"
    else
        mode=single
        echo -ne "$[${total_frags} - ${bad_frags}],"
    fi

    tail -1 $file | grep -oP '^[\d.]+'
done
