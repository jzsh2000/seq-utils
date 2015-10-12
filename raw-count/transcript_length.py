#!/usr/bin/env python

# Date  : 2015-10-09
# Reference: http://seqanswers.com/forums/showthread.php?t=4914

import sys
import HTSeq

gff_file = HTSeq.GFF_Reader(sys.argv[1], end_included=True)

transcripts= {}

for feature in gff_file:
   if feature.type == "exon":
      transcript_id = feature.attr['gene_name']
      if transcript_id == '':
         continue
      if transcript_id not in transcripts:
         transcripts[transcript_id] = list()
      transcripts[transcript_id].append(feature)
      
for transcript_id in sorted(transcripts):      
   transcript_length = 0
   for exon in transcripts[transcript_id]:
      transcript_length += exon.iv.length
   print("%s\t%i\t%i"%(transcript_id, transcript_length, len(transcripts[transcript_id])))
