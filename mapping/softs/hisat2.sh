#!/usr/bin/env bash
set -ue

# output files (in ${outputpath}):
# * ${sampleid}_${lable}.bam
# * ${sampleid}_${lable}.bam.bai
# * ${sampleid}_${lable}.mapping
# * ${sampleid}_${lable}.stat.txt

hisat2=`which hisat2 2>/dev/null || true`
samtools=`which samtools 2>/dev/null || true`
ref=
core=12
fq1=
fq2=
outputpath=`realpath .`
sampleid='sample'
label='hisat2'

# ===== original commands ===================================================
# name=$sampleid_$label
# pathprefix=$outputpath/${sampleid}_${label}
# echo $pathprefix
# 
# $hisat --rna-strandness FR -p $core -x $ref -1 $fq1 -2 $fq2 -S $pathprefix.sam
# $samtools view -@ $core -b -u -S $pathprefix.sam > $pathprefix.bam
# $samtools sort -@ $core $pathprefix.bam -o ${pathprefix}_sort.bam
# mv ${pathprefix}_sort.bam $pathprefix.bam
# $samtools flagstat $pathprefix.bam > $pathprefix.stat.txt
# rm $pathprefix.sam
# 
# $samtools index $pathprefix.bam
# ===========================================================================

usage() {
echo "\
Usage:
    $0 [options] -x <ht2-idx> -1 <m1> [-2 <m2>] [-- [ht2-opt]]
    <ht2-idx>  hisat2 index filename prefix
    <m1>       Files with #1 mates
    <m2>       Files with #2 mates, not used in single-end mode
    <ht2-opt>  Other options passed to hisat2

Options:
    -p/--threads   number of threads (default: 12)
    -i/--id        sampleid (default: sample)
    -l/--label     label (default: hisat2)
    -o/--output    outputpath (default: .)
    --hisat2       hisat2 executable file path
    -s/--samtools  samtools executable file path
    -h/--help      print this usage message
"
}

# -- message color -- #
prompt='\033[1;36m'
error='\033[1;31m'
ncolor='\033[0m'

# ===== read parameters ===== #

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

options=$(getopt \
    -o x:1:2:p:i:l:o:s:h \
    -l threads:,id:,label:,output:,hisat2:,samtools:,help \
    -- $@)
eval set -- "$options"

while true; do
    case "$1" in
        -x)
            ref=$2
            shift 2
            ;;
        -1)
            fq1=$2
            shift 2
            ;;
        -2)
            fq2=$2
            shift 2
            ;;
        -p|--threads)
            core=$2
            shift 2
            ;;
        -i|--id)
            sampleid=$2
            shift 2
            ;;
        -l|--label)
            label=$2
            shift 2
            ;;
        -o|--output)
            outputpath=$2
            shift 2
            ;;
        --hisat2)
            hisat2=$2
            shift 2
            ;;
        -s|--samtools)
            samtools=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift; break
            ;;
        *)
            echo -e "${error}getopt internal error!${ncolor}"
            exit 1
            ;;
    esac
done

# ===== check parameters ===== #

[ $# -gt 0 ] && echo -e "Additional parameters: $@"

# check hisat2
if [ -z "$hisat2" ] || [ ! -x "$hisat2" ]; then
    echo -e "${error}cannot locate hisat2 executable file${ncolor}: $hisat2" >&2
    exit 1
fi
hisat2_version=`$hisat2 --version | head -1 | cut -d' ' -f3`
echo "hisat2 version: $hisat2_version"

# check samtools
if [ -z "$samtools" ] || [ ! -x "$samtools" ]; then
    echo -e "${error}cannot locate samtools executable file${ncolor}: $samtools" >&2
    exit 1
fi
samtools_version=`$samtools 2>&1 | grep -oP '(?<=^Version: )[^ ]+'`
echo "samtools version: $samtools_version"

# check outputpath
if [ ! -d $outputpath ] || [ ! -w $outputpath ]; then
    mkdir -p $outputpath
    if [ ! -d $outputpath ] || [ ! -w $outputpath ]; then
        echo -e "${error}cannot write to outputpath${ncolor}: $outputpath" >&2
        exit 1
    fi
fi

# check sampleid (sequence of [.A-Za-z0-9_])
if ! echo "$sampleid" | grep -qxP '[.\w]+'; then
    echo -e "${error}invalid sampleid${ncolor}: $sampleid" >&2
    exit 1
fi

# check label (sequence of [.A-Za-z0-9_])
if ! echo "$label" | grep -qxP '[.\w]+'; then
    echo -e "${error}invalid label${ncolor}: $label" >&2
    exit 1
fi

# check hisat2-index
if [ ! -f "${ref}".1.ht2 ]; then
    echo -e "${error}cannot load hisat2 index${ncolor}: $ref" >&2
    exit 1
fi

# check fq1 and fq2
if [ -z "$fq1" ] || [ ! -r "$fq1" ]; then
    echo -e "${error}cannot read fastq file 1${ncolor}: $fq1" >&2
    exit 1
fi

if [ -n "$fq2" ] && [ ! -f "$fq2" ]; then
    echo -e "${error}cannot read fastq file 2${ncolor}: $fq2" >&2
    exit 1
fi

[ -z "$fq2" ] && echo "enter single-end mode" || echo "enter paired-end mode"

# check threads
if [ ! $core -ge 1 ]; then
    echo -e "${error}invalid threads number${ncolor}: $core" >&2
    exit 1
fi

# ===== start analysis ===== #
echo -e "${prompt}`date +'%F %T'`\tstart${ncolor}"

name=${sampleid}_${label}
pathprefix=$outputpath/${sampleid}_${label}
echo "writing bam file to: ${pathprefix}.bam"

if [ -z "$fq2" ]; then
    # single-end mode
    $hisat2 $@ -p $core -x $ref -U $fq1 --dta-cufflinks\
        -S $pathprefix.sam \
        2> $pathprefix.mapping
else
    # paired-end mode
    $hisat2 $@ -p $core -x $ref -1 $fq1 -2 $fq2 --dta-cufflinks\
        -S $pathprefix.sam \
        2> $pathprefix.mapping
fi

if [ "${samtools_version:0:1}" == '0' ]; then
    # old version of samtools
    $samtools view -@ $core -b -u -S $pathprefix.sam > $pathprefix.bam
    $samtools sort -@ $core $pathprefix.bam ${pathprefix}_sort
    mv ${pathprefix}_sort.bam $pathprefix.bam
elif [ "${samtools_version:0:3}" == '1.2' ]; then
    $samtools view -@ $core -b -u -S $pathprefix.sam > $pathprefix.bam
    $samtools sort -@ $core $pathprefix.bam -o ${pathprefix}_sort.bam
    mv ${pathprefix}_sort.bam $pathprefix.bam
else
    # newer version of samtools (current 1.3.1)
    $samtools sort -@ $core -o $pathprefix.bam $pathprefix.sam
fi

if [ -f $pathprefix.bam ]; then
    rm $pathprefix.sam

    $samtools index $pathprefix.bam
    $samtools flagstat $pathprefix.bam > $pathprefix.stat.txt
else
    echo -e "${error}cannot find bam file in outputpath${ncolor}" >&2
    exit 1
fi

echo -e "${prompt}`date +'%F %T'`\tend${ncolor}"
