#!/usr/bin/env python3

#converts a fastq to a fasta

import sys

fq=open(sys.argv[1], "r")

for line in fq:
	line=line.strip("\n")
	name=">"+line.strip("@")	
	seq=fq.readline()
	seq=seq.strip("\n")
	a=fq.readline()
	qual=fq.readline()
	print(name, file=sys.stdout)
	print(seq, file=sys.stdout)
	
