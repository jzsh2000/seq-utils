#!/usr/bin/env python3

import pysam
import sys
from argparse import ArgumentParser, FileType

def find_introns(read_iterator):
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
        if 'N' in r.cigarstring:
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

    args = parser.parse_args()
    if not args.bam_file:
        parser.print_help()
        exit(1)

    samfile = pysam.AlignmentFile(sys.argv[1], 'rb')
    res = find_introns(samfile)
    for key,value in res.items():
        print("{}\t{}\t{}\t{}\t{}".format(*key, value))
