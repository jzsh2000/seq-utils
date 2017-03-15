#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

def read_gtf_ss(gtf_ss_file):
    """Return a dictionary {(chr, start, stop): [txs]}
    gtf_ss_file should be a tab-delimited file with following columns:
    * chromosome names
    * start position (not included)
    * end position (not included)
    * strand (+/-)
    * transcript name
    """
    gtf_ss = {}
    while True:
        line = gtf_ss_file.readline()
        if not line: break

        refer, start, end, strand, tx = line.strip().split('\t')
        start, end = int(start), int(end)
        if (refer, start, end, strand) not in gtf_ss:
            gtf_ss[(refer, start, end, strand)] = [tx]
        else:
            gtf_ss[(refer, start, end, strand)].append(tx)

    return gtf_ss

def read_bam_ss(bam_ss_file):
    """Return a dictionary {(chr, start, stop): count}
    gtf_ss_file should be a tab-delimited file with following columns:
    * chromosome names
    * start position (not included)
    * end position (not included)
    * strand (+/-)
    * number of instances
    """
    bam_ss = {}
    while True:
        line = bam_ss_file.readline()
        if not line: break

        refer, start, end, strand, count = line.strip().split('\t')
        start, end, count = int(start), int(end), int(count)
        if end - start > 1:
            bam_ss[(refer, start, end, strand)] = count

    return bam_ss

def check_ss(gtf_ss, bam_ss):
    """Return a concise set of splicing sites
    """
    keys_now = list(bam_ss.keys())
    for key in keys_now:
        if key in gtf_ss:
            continue
        refer, start, end, strand = key
        for start_offset in [-1, 0, 1]:
            for end_offset in [-1, 0, 1]:
                if start_offset == 0 and  end_offset == 0:
                    continue
                key_new = (refer, start+start_offset, end+end_offset, strand)
                if key_new in gtf_ss:
                    bam_ss[key_new] = bam_ss[key] + bam_ss.get(key_new, 0)
                    bam_ss.pop(key)

    return bam_ss

def compare_ss(gtf_ss, bam_ss, count=1):
    """Compare splicing sites in gene annotation and sequence alignment
    print following to standard output:
    - total records in gtf_ss file
    - records in bam_ss file (pass filter)
    - records can be found in gtf_ss file
    - records cannot be found in gtf_ss file
    """
    if count > 1:
        bam_ss_tmp = {x:bam_ss[x] for x in bam_ss if bam_ss[x] >= count}
    else:
        bam_ss_tmp = bam_ss

    print('\t'.join(map(str, [
        len(gtf_ss),
        len(bam_ss_tmp),
        len(set(gtf_ss.keys()).intersection(set(bam_ss_tmp.keys()))),
        len(set(bam_ss_tmp.keys()).difference(set(gtf_ss.keys())))
        ])))

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Compare splice junctions of a BAM file to gene annotation')
    parser.add_argument('gtf_ss',
            nargs=1,
            type=FileType('r'),
            help='input GTF splicing sites')
    parser.add_argument('bam_ss',
            nargs='*',
            type=FileType('r'),
            help='input BAM splicing sites')
    parser.add_argument('-n', '--count',
        dest='count',
        default=1,
        type=int,
        help='only use splicing site with instances >= INT [1]')

    args = parser.parse_args()

    if not args.gtf_ss or not args.bam_ss:
        parser.print_help()
        exit(1)

    gtf_ss = read_gtf_ss(args.gtf_ss[0])

    print('#ss.gtf\tss.bam\tss.known\tss.new')
    for bam_ss_file in args.bam_ss:
        bam_ss = read_bam_ss(bam_ss_file)
        bam_ss = check_ss(gtf_ss, bam_ss)
        compare_ss(gtf_ss, bam_ss, args.count)
