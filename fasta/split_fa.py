#!/usr/bin/env python
# split a fasta file to multiple files

import sys
from Bio import SeqIO

MAX_CHR = 80
records = SeqIO.parse(sys.argv[1], 'fasta')

for seq in records:
    name = seq.name.strip().replace(' ', '_')
    fa_file = open(name, 'w')

    my_seq = []
    for i in range(0, len(seq.seq), MAX_CHR):
        my_seq.append(str(seq.seq[i:i+MAX_CHR]))

    fa_file.write('\n'.join(my_seq))
    fa_file.close()
