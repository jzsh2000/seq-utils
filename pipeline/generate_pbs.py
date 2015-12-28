#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Generate PBS script for RNA-seq mapping and reads counting
# Date: 2015-12-29

import sys
import os
import re

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s DIR\n' % sys.argv[0])
    exit(1)
else:
    script, dname = sys.argv[:2]

(fq1, fq2) = (None, None)
# get fq file
files = os.listdir(dname)
for candidate_file in files:
    if '.fq.gz' in candidate_file:
	if '.1.' in candidate_file:
	    fq1 = os.path.abspath(dname) + '/' + candidate_file
	elif '.2.' in candidate_file:
	    fq2 = os.path.abspath(dname) + '/' + candidate_file

if not fq1 or not fq2:
    sys.stderr.write('Error: cannot locate fastq file in "%s"\n' % dname)
    exit(1)

# get basename of DIR
bname = os.path.split(dname)[-1]
bname = re.sub('[^-A-Za-z0-9_]', '', bname)

if len(bname) < 3:
    bname = 'my_job'

# TODO: estimate phred-score (33 or 64)
print '''
#PBS -S /bin/bash
#PBS -V
#PBS -N %s
#PBS -q batch
#PBS -l nodes=1:ppn=12

symbol='%s'

hisat_cmd='-1 %s -2 %s'

## ================================================== ##

hisat_index='/pfs1/zhanglg/zlg01/Database/hg38/hg38_hisat/hg38_hisat'
gencode_gtf='/pfs1/zhanglg/zlg01/Database/hg38/GENCODE_v23/gencode.annotation.gff3.gz'

cd '%s'
mkdir -p analysis
cd analysis

hisat -p 12 -x $hisat_index ${hisat_cmd} -S ${symbol}.sam &> ${symbol}.hisat.out

samtools view -bSF4 ${symbol}.sam | samtools sort - ${symbol}.sorted
samtools index ${symbol}.sorted.bam
rm ${symbol}.sam

samtools sort -n -@ 12 ${symbol}.sorted.bam ${symbol}.nsorted

htseq-count -i gene_name -s no -f bam ${symbol}.nsorted.bam $gencode_gtf > ${symbol}.count

# vim: set filetype=sh:
''' % (bname, bname, fq1, fq2, os.path.abspath(dname))
