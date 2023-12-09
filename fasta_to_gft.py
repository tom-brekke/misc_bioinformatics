#!/usr/bin/env python
#this script takes a fasta transcriptome (from trinity) as input and another input as the output filename and outputs a .gft file of that filename.
#a gtf file has the list of genes, lengths etc, see the documentation on the web. 

import sys
FileList=sys.argv[1:]
infile=FileList[0] #this is the input trinity fasta file
outfile=FileList[1]  #this is the name of the output file


InFile = open(infile, 'r')
OutFile = open(outfile, 'w')
for line in InFile:
	line=line.strip("\n")
	if line.startswith(">"):	
		fields=line.split(" ")
		gene_ids=fields[0].split(">")[1]
		length=fields[1].split("=")[1]
		comp=gene_ids.split("_seq")[0]
		seq=gene_ids.split("_")[2]
		col1=gene_ids
		col2="Trinity"
		col3="exon"
		col4="1"
		col5=length
		col6="."
		col7="."
		col8="0"
		col9=("%s%s%s%s%s" % ("gene_id \"",comp,"\"; transcript_id \"",gene_ids,"\""))
		newLine=("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (col1, col2, col3, col4, col5, col6, col7, col8, col9))		
		#print newLine		
		OutFile.write(newLine+"\n")
InFile.close()
OutFile.close()		
