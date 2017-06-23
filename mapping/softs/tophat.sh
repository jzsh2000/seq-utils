#!/usr/bin/env bash
set -ue

# output files (in {outputpath})
# .
# ├── accepted_hits.bam
# ├── accepted_hits.bam.bai
# ├── align_summary.txt
# ├── deletions.bed
# ├── insertions.bed
# ├── junctions.bed
# ├── logs/
# ├── prep_reads.info
# └── unmapped.bam

tophat=`which tophat 2>/dev/null || true`
samtools=`which samtools 2>/dev/null || true`
ref=
core=12
fq1=
fq2=
outputpath=`realpath ./tophat_out`

usage() {
echo "\
Usage:
    $0 [options] -x <bt2-idx> -1 <m1> [-2 <m2>] [-- [th-opt]]
    <bt2-idx>  bowtie2 index filename prefix
    <m1>       Files with #1 mates
    <m2>       Files with #2 mates, not used in single-end mode
    <th-opt>   Other options passed to tophat

Options:
    -p/--threads   number of threads (default: 12)
    -o/--output    outputpath (default: ./tophat_out)
    -t/--tophat       tophat executable file path
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
    -o x:1:2:p:o:t:s:h \
    -l threads:,output:,tophat:,samtools:,help \
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
        -t|--tophat)
            tophat=$2
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

# check tophat
if [ -z "$tophat" ] || [ ! -x "$tophat" ]; then
    echo -e "${error}cannot locate tophat executable file${ncolor}: $tophat" >&2
    exit 1
fi
tophat_version=`$tophat --version | cut -d' ' -f2`
echo "tophat version: $tophat_version"

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

# check bowtie2-index
if [ ! -f "${ref}".1.bt2 ]; then
    echo -e "${error}cannot load bowtie2 index${ncolor}: $ref" >&2
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
$tophat $@ -p $core -o $outputpath $ref $fq1 $fq2

if [ -f $outputpath/accepted_hits.bam ]; then
    $samtools index $outputpath/accepted_hits.bam
else
    echo -e "${error}cannot find bam file in outputpath${ncolor}" >&2
    exit 1
fi

echo -e "${prompt}`date +'%F %T'`\tend${ncolor}"

