#!/usr/bin/env Python3


import argparse
import sys
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description="written by TDB\n12-01-2018\nTakes a fasta and calcuates the GC content of each entry as well as GC content of 1st, 2nd, and 3rd codon positions.")                                              
parser.add_argument("--fasta", "-f", type=str, required=True, help="fasta file, typically with ORFs", metavar="*.fasta")
parser.add_argument("--Out", "-o", type=str, required=False, help="Output file name. Omit for stdout")
parser.add_argument("--Out2", "-p", type=str, required=False, help="Output file name for out of phase seqs not to be analyzed. Omit for stdout")
parser.add_argument("--GCxPosition", "-c", default=None, action='store_const', const=True, required=False, help="Include to calculate GC by codon as well as for the entire sequence. Note: Assumes that the sequence is in frame and starts with position 1")
parser.add_argument("--CodonUsage", "-u", default=None, action='store_const', const=True, required=False, help="Include to calculate the codong usage. Note: Assumes that the sequence is in frame and starts with position 1")
args = parser.parse_args()

fasta = args.fasta
Out = args.Out
Out2=args.Out2
#if args.GCxPosition:
#	print("GCxPosition called")
if Out:
	OUT=open(Out,'w')
else:
	OUT=sys.stdout
if Out2:
	OUT2=open(Out2,'w')
else:
	OUT2=sys.stdout

Ns=0 # number of seqs with a non-ACGT character
N=0	#number of sequences that the script analyzed
E=0	#major errors - no stand and no stop  or   start and stop out of sync
EF=0 #frame errors - sequences seems out of frame somehow - has a stop, but no start
usage_dict={"AAA":0, "AAC":0, "AAT":0, "AAG":0, "ACA":0, "ACC":0, "ACT":0, "ACG":0, "ATA":0, "ATC":0, "ATT":0, "ATG":0, "AGA":0, "AGC":0, "AGT":0, "AGG":0, "CAA":0, "CAC":0, "CAT":0, "CAG":0, "CCA":0, "CCC":0, "CCT":0, "CCG":0, "CTA":0, "CTC":0, "CTT":0, "CTG":0, "CGA":0, "CGC":0, "CGT":0, "CGG":0, "TAA":0, "TAC":0, "TAT":0, "TAG":0, "TCA":0, "TCC":0, "TCT":0, "TCG":0, "TTA":0, "TTC":0, "TTT":0, "TTG":0, "TGA":0, "TGC":0, "TGT":0, "TGG":0, "GAA":0, "GAC":0, "GAT":0, "GAG":0, "GCA":0, "GCC":0, "GCT":0, "GCG":0, "GTA":0, "GTC":0, "GTT":0, "GTG":0, "GGA":0, "GGC":0, "GGT":0, "GGG":0,"NNN":0}#any codon with 1 or more "N" will count as NNN
header=False
bases=["A", "C", "G", "T"]


def calc_GC(str, Ns):
	if len(str)==0:
		return("NA", Ns)
	else:	
		str=str.upper()
		GC=str.count("G")+str.count("C")
		AT=str.count("A")+str.count("T")
		N=len(str)-GC-AT
		if N > 0:
			print("Warning: you have some DNA bits with non-ACGT characters. Percent GC is calucated based on the ATGC length only, not including other characters (such as N).", file=sys.stderr)
			#print(str, file=sys.stderr)
			Ns=Ns+1
		#print(str, file=sys.stderr)
		#print(GC, AT, N, len(str), file=sys.stderr)
		perGC=(GC/(GC+AT))*100
		return(perGC, Ns)

def get_in_frame(str, E, EF, name):
	flag=0
	if not str.startswith("ATG"):
		flag=flag+1 
	if not (str.endswith("TAA") or str.endswith("TAG") or str.endswith("TGA")):
		flag=flag+2 
	if len(str) % 3 !=0:
		flag=flag+4
		
	if flag == 0 or flag == 1 or flag == 2 or flag == 6:
		return(str,E,EF)
	elif flag == 3  or flag == 7: #3 should not happen in transcriptomes that are pretty complete. 7 should really never happen at all in this data. 
		if Out2:
			print(">", name, "  flag=", flag, "\n", str, file=OUT2)
		else:
			print("This sequence is really odd - no start or stop - what frame?? Calculate GC by hand! Flag is: ", flag, "\n>", name, "\n", str, sep="", file=sys.stderr)
		E=E+1
		return("",E,EF) #return nothing
	elif flag == 4:
		E=E+1
		if Out2:
			print(">", name, "  flag=", flag, "\n", str, file=OUT2)
		else:
			print("this sequence is really odd - both start or stop are present - but they are out of sync...?? Calculate GC by hand! Flag is: ", flag, "\n", name, "\n", str, sep="", file=sys.stderr)
		return("",E,EF) #return nothing
	elif flag == 5: #ends with stop, but out of frame, add N's to beginning to make up the length 
		EF=EF+1
		if len(str)%3==1:
			return("NN"+str,E,EF)
		elif len(str)%3==2:
			return("N"+str,E,EF)
		elif len(str)%3 ==0:
			print("Warning, this should never occur - the CalcGC script is broken! Exiting now.", file=sys.stderr)
			exit()
			
			
def calc_codon_usage(seq, usage):
	for codonstart in range(0,len(seq),3):
		codon=str(seq[codonstart:codonstart+3])
		if re.search("N", codon):
			usage["NNN"]=usage["NNN"]+1
		else:
			usage[codon]=usage[codon]+1	
	return(usage)	



with open(fasta, 'r') as FASTA:		
	for record in SeqIO.parse(FASTA, 'fasta'):
		usage_dict={"AAA":0, "AAC":0, "AAT":0, "AAG":0, "ACA":0, "ACC":0, "ACT":0, "ACG":0, "ATA":0, "ATC":0, "ATT":0, "ATG":0, "AGA":0, "AGC":0, "AGT":0, "AGG":0, "CAA":0, "CAC":0, "CAT":0, "CAG":0, "CCA":0, "CCC":0, "CCT":0, "CCG":0, "CTA":0, "CTC":0, "CTT":0, "CTG":0, "CGA":0, "CGC":0, "CGT":0, "CGG":0, "TAA":0, "TAC":0, "TAT":0, "TAG":0, "TCA":0, "TCC":0, "TCT":0, "TCG":0, "TTA":0, "TTC":0, "TTT":0, "TTG":0, "TGA":0, "TGC":0, "TGT":0, "TGG":0, "GAA":0, "GAC":0, "GAT":0, "GAG":0, "GCA":0, "GCC":0, "GCT":0, "GCG":0, "GTA":0, "GTC":0, "GTT":0, "GTG":0, "GGA":0, "GGC":0, "GGT":0, "GGG":0,"NNN":0}#any codon with 1 or more "N" will count as NNN
		N=N+1
		#first calc overall GC
		if args.CodonUsage and args.GCxPosition:
			seq,E,EF=get_in_frame(record.seq,E,EF,record.name)
		else:
			seq=record.seq	
		allGC,Ns=calc_GC(seq,Ns)
		#then calc GC from each codon position:
		#print(record.name, record.seq, file=sys.stderr)
		if args.GCxPosition:		
			fir,No=calc_GC(seq[0::3],0)
			sec,No=calc_GC(seq[1::3],0)
			thi,No=calc_GC(seq[2::3],0)
			#here calc the codon usage:
			if args.CodonUsage:
				usage_dict=calc_codon_usage(seq, usage_dict)	
			if header==False: #has the header been added?
				header=True
				print("name", "len_bp", "num_codons", "GC", "GC1", "GC2", "GC3", sep="	", file=OUT, end="	")
				if args.CodonUsage:
					for f in bases:
						for s in bases:
							for t in bases:
								print(f,s,t, sep="", end="	", file=OUT)
					print("NNN", file=OUT, end="")
				print("", file=OUT)					
			print(record.name, len(record.seq), len(record.seq)/3, allGC, fir, sec, thi, sep="	", file=OUT, end="	")
			if args.CodonUsage:
				for f in bases:
						for s in bases:
							for t in bases:
								print(usage_dict[f+s+t], sep="", end="	", file=OUT)
				print(usage_dict["NNN"], file=OUT, end="")
			print("", file=OUT, end="\n")	
		else: #if not args.GCxPosition: #the print statement is above if args.GCxPosition is set
			if header==False:
				header=True
				print("name", "len_bp", "GC", sep="	", file=OUT)
			print(record.name, len(record.seq), allGC, sep="	", file=OUT)
		
if Out:
	OUT.close()	
if Out2:
	OUT2.close()		
	
print(N,"	seqences analyzed", sep="", file=sys.stderr)	
print(E,"	major errors (i.e.: no start and no stop codon, or start and stop out of sync)", sep="", file=sys.stderr)	
print(EF, "	out-of-frame sequences fixed by appending N or NN", file=sys.stderr)
print(Ns, "	sequences with a non-ATG character (may have been introduced when fixing out-of-frame error)", file=sys.stderr)
