#!/usr/bin/env Python3


import argparse
import sys
from Bio import SeqIO
import re
from math import ceil

#print("THIS SCRIPT DOES NOT YET PROPERLY EXTRACT CDS REGIONS")
#quit()


parser = argparse.ArgumentParser(description="written by TDB\n12-01-2018\nTakes a fasta and calcuates the GC content of each entry. Can be used to calculate GC content of 1st, 2nd, and 3rd codon positions using  --GFF. Can calculate GC content for just some sequences using --targets. Can run a window of size X using --Window X and calculate GC for those windows. Mixing and matching some of these functions may work, mixing others may not.")                                              
parser.add_argument("--fasta", "-f", type=str, required=True, help="fasta file, may be simply ORFs or a genome, etc...", metavar="*.fasta")
parser.add_argument("--Out", "-o", type=str, required=False, help="Output file name. Omit for stdout")
parser.add_argument("--Window", "-w", type=int, required=False, help="Used to run a sliding window of the specified size across the sequences. Does not play well with --GFF or --targets.")
parser.add_argument("--targets", "-t", type=str, required=False, help="A list of sequence names to calculate GC content for. If omitted, all sequences will be analyzed.")
parser.add_argument("--GFF", "-g", type=str, required=False, help="Gff3 file with annotations of cds to extract the frame. If omitted, all sequences are assumed to be in frame")
parser.add_argument("--GCxPosition", "-c", default=None, action='store_const', const=True, required=False, help="Include to calculate GC by codon as well as for the entire sequence. Note: Assumes that the sequence is in frame and starts with position 1")
parser.add_argument("--CodonUsage", "-u", default=None, action='store_const', const=True, required=False, help="Include to calculate the codong usage. Note: Assumes that the sequence is in frame and starts with position 1")
parser.add_argument("--override3and7", default=None, action='store_const', const=True, required=False, help="Include to override bailing out on flags 3 and 7. 3 = no start & no stop, 7 = no start, no stop, and not divisible by 3. Normally the GC1, GC2, and GC3 values of these sequences would be NA, using this flag forces the script to treat the seq as in frame as is. Has no effect if GFF and targets file are supplied. Be carefull with this.")
parser.add_argument("--override4", default=None, action='store_const', const=True, required=False, help="Include to override bailing out on flag 4. 4 = start and stop but not divisible by 3. This is probably a mis-annotation and the real stop codon is within the sequence. Normally the GC1, GC2, and GC3 values of these sequences would be NA, using this flag forces the script to search the sequence for the first stop codon in frame from the start codon, and then calculate GC from this fragment. If no stop is found, the entire sequence is returned and GC stats are calculated on the whole thing as if it is in sequence. Has no effect if GFF and targets file are supplied. Probably a good idea to use this if mis-annotatons are common, otherwise leave it out.")
parser.add_argument("--Out2", "-p", type=str, required=False, help="Output file name for out of phase seqs not to be analyzed. Omit for stdout")

args = parser.parse_args()

fasta = args.fasta
Out = args.Out
Out2=args.Out2
#if args.GCxPosition:
#	print("GCxPosition called")
if args.override3and7:
	print("Overridding errors with flags 3 or 7. GC1, GC2, and GC3 will be calculated.", file=sys.stderr)
if args.override4:
	print("Overridding errors with flag 4. GC1, GC2, and GC3 will be calculated.", file=sys.stderr)
if args.Window:
	print("Using a sliding window of size", args.Window, "across all fasta entries", file=sys.stderr)
	W = args.Window

if args.Window and (args.targets or args.GFF or args.CodonUsage):
	print("Some options you have chosen do not play well together. Probably because you are mixing the genome mode (--W) with codon mode (--GFF & --targets & --CodonUsage). Please try again", file=sys.stderr)
	exit()

if Out:
	OUT=open(Out,'w')
else:
	OUT=sys.stdout
if Out2:
	OUT2=open(Out2,'w')
else:
	OUT2=sys.stdout


def calc_codon_usage(seq, usage):
	if len(seq) % 3 ==1:
		seq=seq+"NN"
	elif len(seq) % 3 ==2:
		seq=seq+"N"
	for codonstart in range(0,len(seq),3):
		codon=str(seq[codonstart:codonstart+3])
		if re.search("[BDEFHIJKLMNOPQRSUVWXYZ]", codon):
			usage["NNN"]=usage["NNN"]+1
		else:
			usage[codon]=usage[codon]+1	
	return(usage)	

def calc_GC(str, Ns):
	if len(str)==0:
		return("NA", Ns, 0)
	else:	
		str = str.upper()
		GC = str.count("G") + str.count("C")
		AT = str.count("A") + str.count("T")
		N = len(str) - GC - AT
		if N > 0:
			if Ns==0:
				print("Warning: the input file has at least one sequence with non-ACGT characters. Percent GC is calculated based on the ATGC length only, ignoring other characters such as N and ambiguity codes.", file=sys.stderr)
			#print(str, file=sys.stderr)
			Ns=Ns+1
		#print(str, file=sys.stderr)
		#print(GC, AT, N, len(str), file=sys.stderr)
		if N == len(str):
			perGC = "Undefined"
		else:
			perGC = (GC/(GC+AT))*100
		return(perGC, Ns, N)
		
		
	
	
			
def find_first_stop(str):
	for i in range(0,len(str),3):
		if str[i:i+3].startswith("TAA") or str[i:i+3].startswith("TAG") or str[i:i+3].startswith("TGA"):
			end=i+3
			return(str[0:end])		 
		else:
			return(str)	

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
		if args.override3and7:
			E=E+1
			return(str,E,EF)
		else: # it's unknown what frame these seqs are in - return NA for GC1, GC2, and GC3	
			if Out2:
				print(">", name, "  flag=", flag, "\n", str, file=OUT2)
			else:
				print("This sequence is really odd - no start or stop - what frame?? Calculate GC by hand! Flag is: ", flag, "\n>", name, "\n", str, sep="", file=sys.stderr)
			E=E+1
			return("",E,EF) #return nothing
	
	elif flag == 4:
		E=E+1
		if args.override4:	
			str=find_first_stop(str)
			return(str,E,EF)
		else:
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
			
def get_in_frame_gff(seq, E, EF, name, target_dict_list):
	#print(seq, name, target_dict_list, file=sys.stderr)
	start=int(target_dict_list[0])-1 
	if target_dict_list[1]=="NA": #this happens when a genome is used with multiple GFF gene entires per contig - I force the script to calculate GC for the whole chromosome rather than for each gene.
		stop=len(seq) 
	else:
		stop=int(target_dict_list[1])-1
	frame=int(target_dict_list[2])
	strand=target_dict_list[3]
	length=stop-start
	#print(start,stop,frame,length,length%3,file=sys.stderr)
	if len(seq)<length:
		print("ERROR: GFF sequence length is longer than the fasta sequence length.", "Check that the gff goes with the fasta: ", name, " seq length: ", len(seq),"\nExiting now", file=sys.stderr, sep="")
		exit()
	if frame==0 and start<stop:
		str=seq[start:stop+1:]
	else:
		print("ERROR: sequence is not in frame or is reversed", name, start,stop,frame, file=sys.stderr)
		exit()	
	if strand=="-":
		str=revcom(str)
		#print(str[0:3], str[-3:])	
	flag=0
	#flag explanation: 
	#0/1: has/no start codon
	#0/2: has/no stop codon
	#0/4: is/not multiple of 3 - 
	#0 GREAT! has a start, has a stop, multiple of 3 (in-frame)
	#1 no start, but stop is present and it's a multiple of 3
	#2 no stop, but a start is present
	#3 no start or stop, at least it's a multiple of 3...? Not that that matters in the absense of a start or stop. 
	#4 start and stop, but they can't be in frame b/c length isn't a multiple of 3
	#5 no start and not divisible by 3
	#6 no stop and not divisible by 3
	#7 no start, no stop, not divisible by 3 - this one is a mess. 
	
	if not str.startswith("ATG"):
		flag=flag+1 
	if not (str.endswith("TAA") or str.endswith("TAG") or str.endswith("TGA")):
		flag=flag+2 
	if len(str) % 3 !=0:
		flag=flag+4

	if flag == 0 or flag == 1 or flag == 2  or flag == 6: 
		return(str,E,EF)
	elif flag == 3: #I've mad a new bit for 3: since the annotation assures us that we have the correct frame, but still want it marked that there are a bunch lacking both starts and stops  
		E=E+1
		return(str,E,EF)
	elif flag == 7: #7 should really never happen at all in this data: no start, no stop, not divisible by 3 - this one is a mess.
		if Out2:
			print(">", name, "  flag=", flag, "\n", str, file=OUT2)
		else:
			print("This sequence is really odd - no start or stop - what frame?? Calculate GC by hand! Flag is: ", flag, "\n>", name, "\n", str, sep="", file=sys.stderr)
		E=E+1
		return("",E,EF) #return nothing
	elif flag == 4: #4 start and stop, but they can't be in frame - length isn't a multiple of 3
		E=E+1
		if Out2:
			print(">", name, "  flag=", flag, "\n", str, file=OUT2)
		else:
			print("this sequence is really odd - both start or stop are present - but they are out of sync...?? Calculate GC by hand! Flag is: ", flag, "\n", name, "\n", str, sep="", file=sys.stderr)
		return("",E,EF) #return nothing
	elif flag == 5: #ends with stop, but out of frame, add N's to beginning to make up the length so that it is in frame properly
		EF=EF+1
		if len(str)%3==1:
			return("NN"+str,E,EF)
		elif len(str)%3==2:
			return("N"+str,E,EF)
		elif len(str)%3 ==0:
			print("Warning, this should never occur - the CalcGC script is broken! Exiting now.", file=sys.stderr)
			exit()
			
def get_target_regions(targets):
	target_dict={}
	with open(targets, "r") as TARGETS:
		line = TARGETS.readline().strip()
		bits = line.split("	")
		if length(bits) == 1: #target list is just names of sequences
			target_dict[bits[0]] = 0
			for line in TARGETS:	
				line = line.strip()
				if name not in target_dict:
					target_dict[name] = "all"	
				else:	
					print("ERROR: target list is redundant", file=sys.stderr)

		# elif length(bits) == 3: #target list is a bed file with starts and stops
# 			target_dict[bits[0]] = {bits[1]:bits[2]}
# 			for line in TARGETS:
# 				line = line.strip("\n")
# 				if bits[0] in target_dict:
# 					target_dict[bits[0]][bits[1]] = bits[2] 	
# 				else:
# 					target_dict[bits[0]] = 	{bits[1]:bits[2]}
# 	print(target_dict)						
	return(target_dict)	
	
def parse_GFF(GFFHandle, target_dict):
	with open(GFFHandle, 'r') as GFF:
		for line in GFF:
			if re.search("\tCDS\t", line):
				#print(line)
				line=line.strip("\n")
				fields=line.split("	")
				if fields[0] in target_dict: #only use the ones that are targets
					if target_dict[fields[0]]!=0: #This means there are multiple genes in the GFF file for that contig 
						if args.CodonUsage: #this kills the script - at this point I can't extract multiple genes from each contig.
							print("\n\nThis script does not yet extract genes from contigs with multiple genes - so using '-c' with a genome will not work! \nThe '-c' flag requires one gene/CDS per fasta entry as in a transcriptome. \nExiting now.\n", file=sys.stderr)
							exit()
						else: #force script to calculate GC for the entire contig. 
							start=1
							end="NA" #I don't know the length from the GFF file at this point, but later when this is called, if it's NA, it will be replaced with the entire length of the contig. 
							frame=0
							strand="+"
					else:			
						#print(line)
						start=int(fields[3])
						end=int(fields[4])
						frame=int(fields[7])
						strand=fields[6]			
					target_dict[fields[0]]=start,end,frame,strand
def revcom(seq):
	seq=seq.upper()
	com_dict={"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
	rev=seq[::-1]
	revcom="".join([com_dict[i] for i in rev])
	return(revcom)





Ns=0 # number of seqs with a non-ACGT character
N=0	#number of sequences that the script analyzed
E=0	#major errors - no stand and no stop  or  start and stop out of sync
EF=0 #frame errors - sequences seems out of frame somehow - has a stop, but no start
#usage_dict={"AAA":0, "AAC":0, "AAT":0, "AAG":0, "ACA":0, "ACC":0, "ACT":0, "ACG":0, "ATA":0, "ATC":0, "ATT":0, "ATG":0, "AGA":0, "AGC":0, "AGT":0, "AGG":0, "CAA":0, "CAC":0, "CAT":0, "CAG":0, "CCA":0, "CCC":0, "CCT":0, "CCG":0, "CTA":0, "CTC":0, "CTT":0, "CTG":0, "CGA":0, "CGC":0, "CGT":0, "CGG":0, "TAA":0, "TAC":0, "TAT":0, "TAG":0, "TCA":0, "TCC":0, "TCT":0, "TCG":0, "TTA":0, "TTC":0, "TTT":0, "TTG":0, "TGA":0, "TGC":0, "TGT":0, "TGG":0, "GAA":0, "GAC":0, "GAT":0, "GAG":0, "GCA":0, "GCC":0, "GCT":0, "GCG":0, "GTA":0, "GTC":0, "GTT":0, "GTG":0, "GGA":0, "GGC":0, "GGT":0, "GGG":0,"NNN":0}#any codon with 1 or more "N" will count as NNN
header=False
bases=["A", "C", "G", "T"]

if args.GFF and args.targets:
	target_dict=get_target_regions(args.targets)
	#for k in target_dict:
	#	print(k, target_dict[k])
	parse_GFF(args.GFF, target_dict)
else:
	print("GFF and target file were not both provided, analyzing all sequences.", file=sys.stderr)


with open(fasta, 'r') as FASTA:		
	for record in SeqIO.parse(FASTA, 'fasta'):
		if args.GFF and args.targets: # this bit should skip transcripts that are not in the target list, if a target list was specified
			#print(record.name, target_dict[record.name])
			if record.name not in target_dict:
				#print(record.name, "not in list of targets. Skipping.", file=sys.stderr)
				continue		
		usage_dict={"AAA":0, "AAC":0, "AAT":0, "AAG":0, "ACA":0, "ACC":0, "ACT":0, "ACG":0, "ATA":0, "ATC":0, "ATT":0, "ATG":0, "AGA":0, "AGC":0, "AGT":0, "AGG":0, "CAA":0, "CAC":0, "CAT":0, "CAG":0, "CCA":0, "CCC":0, "CCT":0, "CCG":0, "CTA":0, "CTC":0, "CTT":0, "CTG":0, "CGA":0, "CGC":0, "CGT":0, "CGG":0, "TAA":0, "TAC":0, "TAT":0, "TAG":0, "TCA":0, "TCC":0, "TCT":0, "TCG":0, "TTA":0, "TTC":0, "TTT":0, "TTG":0, "TGA":0, "TGC":0, "TGT":0, "TGG":0, "GAA":0, "GAC":0, "GAT":0, "GAG":0, "GCA":0, "GCC":0, "GCT":0, "GCG":0, "GTA":0, "GTC":0, "GTT":0, "GTG":0, "GGA":0, "GGC":0, "GGT":0, "GGG":0,"NNN":0}#any codon with 1 or more "N" will count as NNN
		N=N+1
		#first calc overall GC
		if args.GFF and args.targets:
			#print(record.name, record.seq, target_dict[record.name], file=sys.stderr)
			seq,E,EF=get_in_frame_gff(record.seq, E, EF, record.name, target_dict[record.name])
			#print(seq,E,EF, file=sys.stderr)
		elif args.CodonUsage and args.GCxPosition:
			seq,E,EF=get_in_frame(record.seq,E,EF,record.name)
		else:
			seq=record.seq
		if args.Window:
			for i in range( ceil( len( record.seq ) / W ) ):
				seq = record.seq[(i*W)+1:(i+1)*W:1] 
				allGC,Ns,numN = calc_GC(seq,Ns)
				if header==False:
					header=True
					print("name", "start", "stop", "num_Ns", "GC", sep="	", file=OUT)
				print(record.name, ( (i * W)+1 ), min(( ( i + 1 ) * W ), len(record.seq)), seq.count("N")+seq.count("n"), allGC, file=OUT)
		else:	
			allGC,Ns,numN = calc_GC(seq,Ns)
			
			#then calc GC from each codon position:
			#print(record.name, record.seq, file=sys.stderr)
			if args.GCxPosition:		
				fir,No,No2 = calc_GC(seq[0::3],1)
				sec,No,No2 = calc_GC(seq[1::3],1)
				thi,No,No2 = calc_GC(seq[2::3],1)
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
				print(record.name, len(seq), len(seq)/3, allGC, fir, sec, thi, sep="	", file=OUT, end="	")
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
					print("name", "len_bp", "NumN","GC", sep="	", file=OUT)
				print(record.name, len(seq), numN, allGC, sep="	", file=OUT)
		
if Out:
	OUT.close()	
if Out2:
	OUT2.close()		
	
print(N,"	seqences analyzed", sep="", file=sys.stderr)	
print(E,"	major errors (i.e.: no start and no stop codon, or start and stop out of sync). If frame was pulled from a gff3, these are probably not an issue.", sep="", file=sys.stderr)	
print(EF, "	out-of-frame sequences fixed by appending N or NN", file=sys.stderr)
print(Ns, "	sequences with a non-ACGT character (may have been introduced when fixing out-of-frame error)", file=sys.stderr)

