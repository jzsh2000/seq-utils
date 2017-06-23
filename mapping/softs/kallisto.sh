#!/usr/bin/env bash
set -ue

# use `kallisto index` to build index

# output files (in ${outputpath}):
# * abundance.h5
# * abundance.tsv (main output)
# * run_info.json

kallisto=`which kallisto 2>/dev/null || true`
ref=
core=1
fq1=
fq2=
outputpath=`realpath ./kallisto_out`

usage() {
echo "\
Usage:
    $0 [options] -x <index> -1 <m1> [-2 <m2>] [-- [options]]
    <index>    Filename for the kallisto index
    <m1>       Fastq files with #1 mates
    <m2>       Fastq files with #2 mates, not used in single-end mode
    <options>  Other options passed to BWA

Options:
    -p/--threads   number of threads (default: 1)
    -o/--output    outputpath (default: .)
    -k/--kallisto  kallisto executable file path
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
    -o x:1:2:p:o:h \
    -l threads:,output:,help \
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
        -o|--output)
            outputpath=$2
            shift 2
            ;;
        -k|--kallisto)
            kallisto=$2
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

# check outputpath
if [ ! -d $outputpath ] || [ ! -w $outputpath ]; then
    mkdir -p $outputpath
    if [ ! -d $outputpath ] || [ ! -w $outputpath ]; then
        echo -e "${error}cannot write to outputpath${ncolor}: $outputpath" >&2
        exit 1
    fi
fi

# check kallisto index
if [ ! -r "${ref}" ]; then
    echo -e "${error}cannot load kallisto index${ncolor}: $ref" >&2
    exit 1
fi

# check kallisto executable file
if [ ! -x "${kallisto}" ]; then
    echo -e "${error}cannot locate kallisto executable file${ncolor}: $kallisto" >&2
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
kallisto quant $@ -t ${core} -i ${ref} -o ${outputpath} $fq1 $fq2
echo -e "${prompt}`date +'%F %T'`\tend${ncolor}"
