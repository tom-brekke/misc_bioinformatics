#!/usr/bin/env python
#converts a tab-delimted "fasta" back to a regular fasta" i.e., does the opposite of Convert_Fasta_to_Tab.pl
#takes two inputs, an infile and an outfile
#for each line of infile
#replace the second tab with a \n
#write to outfile name

import sys
FileList=sys.argv[1:]
infile=FileList[0] 
outfile=FileList[1]

IN=open(infile, "r")
OUT=open(outfile, "w")
for line in IN:
	line=line.strip("\n")
	fields=line.split("	")
	name=fields[0]
	seq=fields[2]
	lineA=("%s%s%s%s" % (">",name," ",fields[1]))	
	OUT.write(lineA+"\n"+seq+"\n")
	
IN.close()
OUT.close()
