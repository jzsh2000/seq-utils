#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Generate PBS script for RNA-seq mapping and reads counting
# Date: 2015-12-29

import sys
import os
import re
import gzip

if len(sys.argv) < 2:
    # should be one sample per directory
    sys.stderr.write('Usage: %s DIR\n' % sys.argv[0])
    exit(1)
else:
    script, dname = sys.argv[:2]

inc_file_path = os.path.dirname(os.path.abspath(script)) + '/inc/path'
inc_file = open(inc_file_path, 'r')
if not inc_file:
    sys.stderr.write('cannot locate file %s\n' % inc_file_path)

# read in config file
path_dict = {}
while True:
    line = inc_file.readline().strip()
    if not line: break

    line = line.split(':')
    path_dict[line[0]] = line[1]

inc_file.close()

hisat_index = path_dict['hisat_index']
count_script = path_dict['count_script']
count_script_se = path_dict['count_script_se']

def get_fastq(dname):
    '''get fq.gz files'''
    (fq0, fq1, fq2) = (None, None, None)
    files = [ x for x in os.listdir(dname) if '.fq.gz' in x ]
    for candidate_file in files:
        fq_abspath = os.path.abspath(dname) + '/' + candidate_file
        if re.search('[-._]1\.', candidate_file):
            fq1 = fq_abspath
        elif re.search('[-._]2\.', candidate_file):
            fq2 = fq_abspath
        else:
            fq0 = fq_abspath

    if fq0:
        return (fq0,)

    if fq1 and fq2:
        return (fq1, fq2)
    elif fq1:
        fq0 = fq1
        return (fq1,)
    elif fq2:
        fq0 = fq2
        return (fq2,)

    sys.stderr.write('Error: cannot locate fastq file in "%s"\n' % dname)
    exit(1)

def get_phred_score(fq):
    fq_file = gzip.open(fq, 'rb')
    # read in some reads and check the quality score
    quality_list = fq_file.readlines(100000)[3::4]
    fq_file.close()

    if re.search('[O-h]', ''.join(quality_list)):
        return 64
    else:
        return 33

def get_job_symbol(dname):
    '''get basename of DIR'''
    bname = os.path.split(dname)[-1]
    bname = re.sub('[^-A-Za-z0-9_.]', '', bname)

    if len(bname) < 3:
        bname = 'my_job'

    return(bname)

def write_script(bname, fqs, phred_score):
    '''print pbs script to stdout'''
    if(len(fqs) == 1):
        mode = 'se'
    else:
        mode = 'pe'

# TODO: estimate phred-score (33 or 64)
    print '''#PBS -S /bin/bash
#PBS -V
#PBS -N %s
#PBS -q batch
#PBS -l nodes=1:ppn=12

symbol='%s'
    ''' % (bname, bname)

    if mode == 'pe':
        print "hisat_cmd='-1 %s -2 %s'" % (fqs[0], fqs[1])
    else:
        print "hisat_cmd='-U %s'" % fqs[0]

    print '''
## ================================================== ##

hisat_index='%s'
count_script='%s'
count_script_se='%s'
''' % (hisat_index, count_script, count_script_se)

    if(phred_score == 64):
        print "phred_option='--phred64'"
    else:
        print "phred_option='--phred33'"

    print '''
cd '%s'
mkdir -p fastqc bam count
''' % os.path.abspath(dname)

    if mode == 'pe':
        print "fastqc -t 2 -o 'fastqc/' %s %s" % (fqs[0], fqs[1])
    else:
        print "fastqc -o 'fastqc/' %s" % (fqs[0])

    print '''
cd bam
hisat -p 12 -x $hisat_index $phred_option ${hisat_cmd} -S ${symbol}.sam &> ${symbol}.hisat.out

samtools view -bSF4 ${symbol}.sam | samtools sort -@ 12 - ${symbol}.sorted
samtools index ${symbol}.sorted.bam
samtools sort -n -@ 12 ${symbol}.sorted.bam ${symbol}.nsorted
rm ${symbol}.sam
cd ..
'''
    if mode == 'pe':
        print "$count_script bam/${symbol}.sorted.bam > count/${symbol}.count"
    else:
        print "$count_script_se bam/${symbol}.sorted.bam > count/${symbol}.count"

    print '\n# vim: set filetype=sh:'
    pass

fqs = get_fastq(dname)
phred_score = get_phred_score(fqs[0])
bname = get_job_symbol(dname)
write_script(dname, fqs, phred_score)
