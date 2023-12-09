#!/usr/bin/env python3	
#created by TDB 1-27-2015
#a script to make a gtf out of a transcriptome

import sys
import re

tran_file=sys.argv[1]
gtf_file=sys.argv[2]

tran=open(tran_file, "r")
gtf=open(gtf_file, "w")

#print("##gff-verson 3", file=gtf)
for line in tran:
	if line.startswith(">"):
		fields=re.search(">(comp\d+_c\d+_seq\d+)\s+len=(\d+)", line)
		seqname=fields.group(1)
		source="Trinity"
		feature="isoform"
		start=1
		end=fields.group(2)
		score="."
		strand="."
		frame="."
		attr="transcript_id \""+seqname+"\";gene_id \""+seqname.split("_seq")[0]+"\""
		print(seqname, source, feature, start, end, score, strand, frame, attr, sep="\t", file=gtf)
