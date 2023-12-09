#!/usr/bin/env Python3
import sys
from Bio import SeqIO


usage='''
written by TDB
19-07-2017
a script to pull out genomic regions of a reference, based on the alignments in a sam file. 

sys.argv inputs are:
	#[1] the genome fasta:		GCA_002204375.1_MunDraft-v1.0_genomic.fna
	#[2] the sam file: 		GBS_Lucigen_alignment.sam
	#[3] the vcf:			batch_1_reformatted.vcf
'''
if len(sys.argv)<4:
	print(usage, file=sys.stderr)
	exit()

ODDITIES=open("ODD.fasta", "w")
WELLBEHAVED=open("GOOD.fasta", "w")
NEWVCF=open("batch_1_reformatted_with_Lucigen.vcf", "w")


recognitionSite="CAGCT" #the cut site that GBS reads should start with.


#######definitions
def countMatches(s, t):
	counter=0
	for i in range(0,len(s)):
		if s[i] == t[i]:
			counter+=1
	return(counter)

def getMismatches(s, t):
	mismatch_list=[]
	for i in range(0,len(s)):
		if s[i] != t[i]:
			mismatch_list.append(tuple((i,s[i],t[i]))) #this tuple is the position, then the basecalls of s and t
	return(mismatch_list)
		
		
def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N':'N', 'n':'n'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return(''.join(letters))
def revcom(s):
    return(complement(s[::-1]))
#############



#############
#make a dictionary where keys are genome contig names and values are a list of locations
dict={}
print("\nmaking dictionary", file=sys.stderr)
with open(sys.argv[2], "r") as SAM:
	for line in SAM:
		line=line.strip("\n")
		if not line.startswith("@"):
			fields=line.split("	")
			QNAME=fields[0]
			FLAG=int(fields[1])
			RNAME=fields[2]
			POS=int(fields[3])
#			MAPQ=fields[4]
			CIGAR=fields[5]
#			RNEXT=fields[6]
#			PNEXT=fields[7]
#			TLEN=fields[8]
			GSEQ=fields[9]
#			QUAL=fields[10]
			if RNAME not in dict:
				dict[RNAME]=[(QNAME, POS, FLAG, CIGAR, GSEQ)]	
			else:
				dict[RNAME].append((QNAME,POS, FLAG, CIGAR, GSEQ))
print("  -finished making dictionary\n", file=sys.stderr)
#############


#############
#make the dictionary of the VCF:
print("reading VCF", file=sys.stderr)
vcf_dict={}
with open(sys.argv[3], "r") as VCF:
	for line in VCF:
		line=line.strip("\n")
		if not line.startswith("#"):
			fields=line.split("	")
			TAG=fields[0]
			if TAG not in vcf_dict:		
				vcf_dict[TAG]=[tuple(fields)]
			else:
				vcf_dict[TAG].append(tuple(fields))
print("  - finished reading VCF", file=sys.stderr)
#############




#############
#extract the appropriate genomic regions
print("extracting genomic regions", file=sys.stderr)
#second, read the genome record by record and get the bits needed:
with open(sys.argv[1], "r") as GENOME:
	for record in SeqIO.parse(GENOME, "fasta"):
#		check if record.id is in the dictionary to extract and if so, do the extraction
		if record.id in dict:
			for ENTRY in dict[record.id]:
				QNAME=ENTRY[0]
				START=ENTRY[1]-1
				END=START+92
				FLAG=ENTRY[2]
				CIGAR=ENTRY[3]
				GSEQ=ENTRY[4] 				#GBS-SEQ
				LSEQ=record.seq[START:END]	#Lucida genome seq
				if FLAG == 16:
					LSEQ=revcom(LSEQ)
					GSEQ=revcom(GSEQ)	
				if CIGAR=="92M" and countMatches(LSEQ.upper()[0:len(recognitionSite)], recognitionSite) >= len(recognitionSite)-1: # the '>=len()-1' bit allows 1 mismatch		
					misMatches=getMismatches(GSEQ, LSEQ.upper())
					#print(misMatches)
					if len(misMatches) >0:
						#get record from vcf:
						if QNAME in vcf_dict:
							print(vcf_dict[QNAME])
							print(misMatches, "\n")
							#check if SNP is in vcf_dict by stepping through SNPs and checking each one.
							for mis in misMatches:
								flag=0 #this keeps track of whether the mismatch is in the vcf
								for record in vcf_dict[QNAME]:
									if mis[0] == record[1]:#mis[0] is the first entry of the tuple from getMismatches which is the position. record[1] is the position bit of the vcf record
										flag=1
										#add a record for the Lucinda
										addRecord		
								if flag==0:
								
									#add a whole new record for everyone.		
						#else: #this occurs when there is a SNP in a tag that didn't used to have any snps - I need to make a new record 		
								#make a new record
					print(QNAME,"	", LSEQ,sep="", file=WELLBEHAVED)
					
				else:
					print(QNAME,"	", LSEQ,sep="", file=ODDITIES)
					
print("  -finished extracting genomic regions\n", file=sys.stderr)

			
ODDITIES.close()
WELLBEHAVED.close()					
NEWVCF.close()			