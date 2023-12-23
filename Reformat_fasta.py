#!/usr/bin/env Python3

#this script was written by TDB on 07-09-2016
#it takes in a fasta file and outputs another fasta file with the sequence all on one line.

import sys
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description="written by TDB\n07-09-2017\nReformats a *.fasta file so that the sequence is all on one line. Also can generate a *.bed or a *.gtf from any *.fasta input with or without outputting the 1-line fasta.")                                              
parser.add_argument("-F", "--FASTA", type=str, required=True, help="A fasta file.", metavar="*.fasta")
parser.add_argument("-o", "--out", type=str, required=False, help="Only used if -f1 is called, otherwise output will go to stdout. The output 1 line fasta name. Default stdout", default=sys.stdout, metavar="*.fasta")
parser.add_argument("-f1", "--OneLineFasta", type=bool, required=False, help="Make a 1 line fasta file? True or False. Will go to stdout unless -o is set. Default: %(default)s", default=False, metavar="True")
parser.add_argument("-b", "--BED", type=bool, required=False, help="Make a bed file? True or False. Default: %(default)s", default=False, metavar="True")
parser.add_argument("-g", "--GTF", type=bool, required=False, help="Make a gtf file? True or False. Default: %(default)s", default=False, metavar="False")

args = parser.parse_args()

FASTA = args.FASTA
OUT = args.out
makeBED = args.BED
makeGTF = args.GTF
makeOLFasta = args.OneLineFasta


fileEnd=FASTA.split(".")[-1]

#prepare the files to write out to:
if makeOLFasta:
	if OUT != sys.stdout:
		OUT=open(OUT, 'w')

if makeBED:
	BEDname=FASTA.replace(fileEnd, "bed")
	BEDOUT=open(BEDname, 'w')

if makeGTF:
	GTFname=FASTA.replace(fileEnd, "gtf")
	GTFOUT=open(GTFname, 'w')
		




#definitions for making various things:
def makeOLFastafile(record):
	print(record.id, file=OUT)
	print(record.seq, file=OUT)

def makeBEDfile(record):
	col1=record.id	#The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
	col2="0"	#The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
	col3=len(record.seq)	#chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
	print(col1, col2, col3, sep="	",  end="\n", file=BEDOUT)

def makeGTFfile(record): # a pretty simple GTF file - only name, start, and stop. Added this code mainly as a placeholder in case of future improvements:
	col1=record.id	#seqname
	col2="."	#source: eg. Ensembl
	col3="."	#feature: eg. Gene, Variation
	col4="1"	#start position, first at 1 not 0
	col5=len(record.seq)	#end position first at 1 not 0
	col6="."	#score: floating point
	col7="."	#strand: +,-
	col8="."	#frame: 0,1,2
	col9="."	#attribute semicolon separated list of tag-value pairs
	print(col1, col2, col3, col4, col5, col6, col7, col8, col9, sep="	", end="\n", file=GTFOUT)





#body of the code:
for record in SeqIO.parse(FASTA, "fasta"):
		if makeOLFasta:
			makeOLFastafile(record)
		if makeBED:
			makeBEDfile(record)
		if makeGTF:
			makeGTFfile(record)
		
		
						 
#close the files.
if makeOLFasta:
	OUT.close() 
if makeBED:
	BEDOUT.close()
if makeGTF:
	GTFOUT.close()
