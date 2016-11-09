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
        if 'N' in r.cigarstring and r.mapping_quality>=mapq_threshold:
            last_read_pos = False
            for read_loc, genome_loc in r.get_aligned_pairs():
                refer = r.reference_name
                strand = '-' if r.is_reverse else '+'
                if read_loc is None and last_read_pos:
                    start = genome_loc
                elif read_loc and last_read_pos is None:
                    stop = genome_loc  # we are right exclusive ,so this is correct
                    res[(refer, start, stop, strand)] += 1
                    del start
                    del stop
                last_read_pos = read_loc
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
