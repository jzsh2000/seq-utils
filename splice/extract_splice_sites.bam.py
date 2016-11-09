#!/usr/bin/env python3

import pysam
import sys
from argparse import ArgumentParser, FileType

def find_introns(read_iterator, mapq_threshold):
    """Return a dictionary {(chr, start, stop): count}
    Listing the intronic sites in the reads (identified by 'N' in the cigar strings),
    and their support ( = number of reads ).

    read_iterator can be the result of a .fetch(...) call.
    Or it can be a generator filtering such reads. Example
    find_introns((read for read in samfile.fetch(...) if read.is_reverse)
    """
    
    import collections
    res = collections.Counter()
    for r in read_iterator:
        if r.mapping_quality >= mapq_threshold and 'N' in r.cigarstring:
            refer = r.reference_name
            strand = '-' if r.is_reverse else '+'
            exons = r.get_blocks()
            exons.sort()
            for i in range(1, len(exons)):
                res[(refer, exons[i-1][1], exons[i][0]+1, strand)] += 1
    return res

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract splice junctions from a BAM file')
    parser.add_argument('bam_file',
        nargs='?',
        type=FileType('r'),
        help='input BAM file')
    parser.add_argument('-q', '--mapq',
        dest='mapq',
        default=30,
        type=int,
        help='only include reads with mapping quality >= INT [30]')

    args = parser.parse_args()
    if not args.bam_file:
        parser.print_help()
        exit(1)

    samfile = pysam.AlignmentFile(args.bam_file, 'rb')
    res = find_introns(samfile, args.mapq)
    for key,value in res.items():
        print("{1}\t{2}\t{3}\t{4}\t{0}".format(value, *key))
