#!/usr/bin/env Python3


import argparse
import sys
from Bio import SeqIO
#import gffutils
#from BCBio.GFF import GFFParser,GFFExaminer
import re
from math import ceil
import pprint


parser = argparse.ArgumentParser(description="written by TDB\n02-09-2019\nTakes a fasta and a gff. Extracts every CDS, stitches them together into mRNAs and calcuates the GC content of each one. Can be used to calculate GC content of 1st, 2nd, and 3rd codon positions as well.")                                              
parser.add_argument("--fasta", "-f", type=str, required=True, help="fasta file of a genome", metavar="*.fasta")
parser.add_argument("--GFF", "-g", type=str, required=False, help="Gff3 file with annotations of cds to extract the frame. If omitted, all sequences are assumed to be in frame")
parser.add_argument("--Output_type", "-t", type=str, required=True, help="Type of output: can be either 'fasta' or 'GC_table'")
parser.add_argument("--GCxPosition", "-c", default=None, action='store_const', const=True, required=False, help="Include to calculate GC by codon as well as for the entire sequence. Note: Assumes that the sequence is in frame and starts with position 1")
parser.add_argument("--Out", "-o", type=str, required=False, help="Output file name. Omit for stdout")

args = parser.parse_args() 

if args.Out:
	OUT=open(args.Out,'w')
else:
	OUT=sys.stdout



def add_exon_or_cds_to_dict(bits, chr):
	#print("cds", bits, file=sys.stderr)
	source     = bits[1]
	type       = bits[2] 
	start      = int(bits[3])
	end        = int(bits[4])
	score      = bits[5]
	strand     = bits[6]
	frame      = bits[7]
	info       = bits[8]	
	Name       = parse_info(info, "Gene")
	ID         = parse_info(info, "ID")	
	parent_ID  = parse_info(info, "parent")
	try:
		parent_type = parent_type_dict[parent_ID]
	except:
		print("Error 1. This error may be thrown because a gene crossed the boundary of an ECR - grep the appropriate ECR gff for this:\n	", parent_ID, "\nand consider revising the exact breakpoints be a couple bases for that species if it seems reasonable. \nExiting the subroutine now.", file=sys.stderr, sep="")
		#specifically, the error happens when an exon is assigned to an ECR while it's parent gene/mRNA was assigned to the space between two ECRs (or another ECR)
		#easiest way to fix it is to shift the ECR breakpoints so that they don't span a gene. May not be possible in all cases if the gene spans the entire breakpoint and exists in two separate ECRs..
		return("nope")
	#check to make sure that the parent is indeed an mRNA - it may be a pseudogene, in which case I want to just skip it.
	if parent_type =="mRNA": #if the parent is a gene, go ahead 
		if parent_ID in parent_tree_dict: #this will end up skipping any cds/exon that has as it's parent a gene in stead of an mRNA - I've only seen this is cases where there are V_gene_segements which are not super common. 
			grandparent_ID = parent_tree_dict[parent_ID]#this comes from looking up the tree to the parent RNA... 
			#print(chr, grandparent_ID, parent_type, parent_ID, file=sys.stderr)
			if grandparent_ID not in gff_dict[chr]["gene"]: #this means that the gene spanned a breakpoint region. The 'gene' part of the record is almost certainly under the key for the original Chromosome. Need to search it out and add it in to the appropriate ECR - this may end up adding each gene to multiptle ECRs depending on whether it crossed the entire breakpoint and into the ECR nextdoor, or just reached into the breakpoint. Refining breakpoints to avoid splitting genes may be a good approach to get around this. 
				for tmpChr in gff_dict:
					if grandparent_ID in gff_dict[tmpChr]["gene"]:
						#print(tmpChr, chr, grandparent_ID, file=sys.stderr)
						gff_dict[chr]["gene"][grandparent_ID] = gff_dict[tmpChr]["gene"][grandparent_ID] #this moves over the entry. But there can still be an error when
					else:
						pass
						#print(chr, tmpChr, grandparent_ID, parent_ID, ID)
						#print("Error - a gene that seems to cross a breakpoint and can't be found in the gff dict. This may be a tricky error to sort out. The 'gene' bit of the gff was never added to the gff_dict or was added under the original chr rather than the ecr. \nGood luck. \nExiting now.", file=sys.stderr)
						#exit()	
			#print(chr, "gene", grandparent_ID, parent_type, parent_ID)
			if parent_ID not in gff_dict[chr]["gene"][grandparent_ID][parent_type]:
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID] = {"Name":Name, "type":type, "source":source, "start":start, "end":end, "score":score, "strand":strand, "frame":frame, "info":info, "subfeatures":{}}
			if type not in gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"]:#make sure that both "exon " and "cds" are headings under mRNA
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type] = {}
			if ID not in gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type]:
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type][ID] = {"Name":Name, "type":type, "source":source, "start":[start], "end":[end], "score":[score], "strand":[strand], "frame":[frame], "info":info}	
			else:# if the ID is already in there, it's probably a CDS (all have the same ID) and the starts etc just need to be updated
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type][ID]["start"].append(start)
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type][ID]["end"].append(end)
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type][ID]["strand"].append(strand)
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type][ID]["frame"].append(frame)
				gff_dict[chr]["gene"][grandparent_ID][parent_type][parent_ID]["subfeatures"][type][ID]["score"].append(score)


def add_gene_or_region_to_dict(bits, chr):
	#print("gene", bits, file=sys.stderr)
	source     = bits[1]
	type       = bits[2]
	start      = int(bits[3])
	end        = int(bits[4])
	score      = bits[5]
	strand     = bits[6]
	frame      = bits[7]
	info       = bits[8]
	ID         = parse_info(info, "ID")	
	Name       = parse_info(info, "Name")
	if ID not in gff_dict[chr][type]:
		gff_dict[chr][type][ID] = {"Name":Name, "type":type, "source":source, "start":start, "end":end, "score":score, "strand":strand, "frame":frame, "info":info, "mRNA":{}, "V_gene_segment":{}, "C_gene_segment":{}}		



def add_mRNA_to_dict(bits, chr):
	#print("mRNA", bits, file=sys.stderr)
	source     = bits[1]
	type       = bits[2] 
	start      = int(bits[3])
	end        = int(bits[4])
	score      = bits[5]
	strand     = bits[6]
	frame      = bits[7]
	info       = bits[8]
	ID         = parse_info(info, "ID")
	Name       = parse_info(info, "Gene")
	parent_ID  = parse_info(info, "parent")
	if ID not in parent_tree_dict: # set up the mRNA_ID:GeneID dictionary for later lookup
		parent_tree_dict[ID] = parent_ID
	if parent_ID not in gff_dict[chr]["gene"]: #if the parent is in the pseudogenes, that means there is no mRNA and the next bit should be skipped. If it does not exist there, then it is an mRNA and we're good to go. 
		gff_dict[chr]["gene"][parent_ID] = {"Name":parent_ID, "type":"NA", "source":source, "start":-1, "end":-1, "score":0, "strand":"NA", "frame":-1, "info":"Unknown", "mRNA":{}, "V_gene_segment":{}, "C_gene_segment":{}} #add it in
	if type not in gff_dict[chr]["gene"][parent_ID]:
		gff_dict[chr]["gene"][parent_ID][type] = {}
	if ID not in gff_dict[chr]["gene"][parent_ID][type]:
		gff_dict[chr]["gene"][parent_ID][type][ID] = {"Name":Name, "type":type, "source":source, "start":start, "end":end, "score":score, "strand":strand, "frame":frame, "info":info, "subfeatures":{}}


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
			Ns=Ns+1
		if N == len(str):
			perGC = "Undefined"
		else:
			perGC = ( GC / ( GC + AT )) * 100
		return(perGC, Ns, N)



def extract_cds(output_type):
	Ns = 0
	header=True
	if output_type == "fasta":
		print("Making fasta:", args.OUT, file=sys.stderr)
	elif output_type == "GC_table":
		print("Making GC table", file=sys.stderr)
	for chr in gff_dict:
		for gene in gff_dict[chr]["gene"]:
			#print(chr, gene, file=sys.stderr)
			if "mRNA" in gff_dict[chr]["gene"][gene]:
				for mRNA in gff_dict[chr]["gene"][gene]["mRNA"]:
					if "CDS" in gff_dict[chr]["gene"][gene]["mRNA"][mRNA]["subfeatures"]:
						for CDS in gff_dict[chr]["gene"][gene]["mRNA"][mRNA]["subfeatures"]["CDS"]:
							starts  = gff_dict[chr]["gene"][gene]["mRNA"][mRNA]["subfeatures"]["CDS"][CDS]["start"]
							ends    = gff_dict[chr]["gene"][gene]["mRNA"][mRNA]["subfeatures"]["CDS"][CDS]["end"]
							strands = gff_dict[chr]["gene"][gene]["mRNA"][mRNA]["subfeatures"]["CDS"][CDS]["strand"]
							seq = extract_RNA_seq(chr, starts, ends, strands)
							if output_type == "fasta":
								print("Outputting fasta record for chr:", chr, ", gene:", gene, ", mRNA:", mRNA, ", CDS:", CDS, sep="", end="\r", file=sys.stderr)
								print(">",gff_dict[chr]["gene"][gene]["Name"], "CHR"+chr, "Start"+starts[0], "GeneID"+gene, "mRNA"+mRNA, "CDS:"+CDS,"\n", seq, file=OUT)# to print as a fasta
							elif output_type == "GC_table":
								print("Calculating GC for chr:", chr, ", gene:", gene, ", mRNA:", mRNA, ", CDS:", CDS, sep="", end="\r", file=sys.stderr)
								allGC,Ns,numN = calc_GC(seq, Ns)
								flag = check_assembly(seq)
								if args.GCxPosition:
									fir,No,No2 = calc_GC(seq[0::3],1)
									sec,No,No2 = calc_GC(seq[1::3],1)
									thi,No,No2 = calc_GC(seq[2::3],1)
									if header: #has the header been added?
										header = False
										print("name", "chr", "start", "GeneID", "mRNA", "CDS", "len_bp", "num_codons", "flag", "NumN", "GC", "GC1", "GC2", "GC3", sep="	", file=OUT, end="\n")
									print(gff_dict[chr]["gene"][gene]["Name"], chr, starts[0], gene, mRNA, CDS, len(seq), round(len(seq)/3, 0), flag, numN, allGC, fir, sec, thi, sep="	", file=OUT, end="\n")
								else:
									if header:
										header = False
										print("name", "chr", "start", "GeneID", "mRNA", "CDS", "len_bp", "flag", "NumN","GC", sep="	", file=OUT, end="\n")
									print(gff_dict[chr]["gene"][gene]["Name"], chr, starts[0], gene, mRNA, CDS, len(seq), flag, numN, allGC, sep="	", file=OUT)
	print("", file=sys.stderr)


def extract_RNA_seq(chr, starts, ends, strands):	
	seq = ""
	for s,e,d in zip(starts, ends, strands):
		if d == "+":
			seq = seq + record_dict[chr].seq[s-1:e]			
		elif d == "-":
			seq = seq + revcom(record_dict[chr].seq[s-1:e])
	#check for being in frame:
	return(seq)	


def check_assembly(str):
	flag=0
	#flag explanation: 
	#0/1: has/no start codon
	#0/2: has/no stop codon
	#0/4: is/not multiple of 3 - 
	#so, in sum:
	#0 GREAT! has a start, has a stop, multiple of 3 (in-frame)
	#1 no start, but stop is present and it's a multiple of 3
	#2 no stop, but a start is present
	#3 no start or stop, at least it's a multiple of 3...? Not that that matters in the absense of a start or stop. 
	#4 start and stop, but they can't be in frame b/c length isn't a multiple of 3
	#5 no start and not divisible by 3
	#6 no stop and not divisible by 3
	#7 no start, no stop, not divisible by 3 - this one is a mess. 
	if not str.startswith("ATG"):
		flag = flag + 1 
	if not (str.endswith("TAA") or str.endswith("TAG") or str.endswith("TGA")):
		flag = flag + 2 
	if len(str) % 3 != 0:
		flag = flag + 4
	return(flag)

def parse_fasta(fasta_in_file):
	print("Parsing genome:", fasta_in_file, file=sys.stderr)
	with open(fasta_in_file, "r") as FASTA:
		for record in SeqIO.parse(FASTA, 'fasta'):
			print("	Reading chr", record.id, "                           ", file=sys.stderr, end="\r")
			if record.id not in record_dict:
				record_dict[record.id] = record
			else:
				print(record.id, "already in record_dict. This is odd. Exiting now.", file = sys.stderr)
				exit()


def parse_gff(gff_in_file):
	print("Parsing GFF file:", gff_in_file, file=sys.stderr)
	with open(gff_in_file, 'r') as FH:
		for line in FH:
			if not line.startswith("#"):
				bits = line.strip("\n").split("	")
				chr        = bits[0]
				type       = bits[2]
				info       = bits[8]
				ID         = parse_info(info, "ID")
				if ID not in parent_type_dict:
					parent_type_dict[ID] = type
				if chr not in record_dict:# only add in things that have a fasta record associated, kill script otherwise. Perhaps eventually, I could have them add a new seq_record object to the dict, but not for now. 
					if not re.search("ECR_breakpoint_region", chr): 
						if chr not in alias_dict and type == "region":#this deals with the peromyscus aliases
							alias_dict[chr] = parse_info(bits[8], "Alias")
						chr = alias_dict[chr]
						#print(line, "\n, chr is", chr, file=sys.stderr)
						if not chr is None and chr not in record_dict:
							print("Error, a feature in the gff does not have a corresponding entry in the fasta file!\nExiting now.", file=sys.stderr)
							#exit()
					else:#this means that the gene is in a chromosomal breakpoint region
						continue	
			#this means there is a fasta sequence record that I will be able to use. 
				if type in ["region", "gene", "pseudogene", "lnc_RNA", "tRNA", "rRNA", "pseudogenic_transcript", "ncRNA_gene"]: # these are higher-level features, annotate the top-level bit				
					add_gene_or_region_to_dict(bits, chr)	
				elif type in ["mRNA", 'transcript', "V_gene_segment", "C_gene_segment", "snRNA", "miRNA", "snoRNA", "ncRNA", "scRNA", "J_gene_segment", "D_gene_segment", "gene_segment", "five_prime_UTR", "three_prime_UTR"]: # this is a mid-level feature, add in the gene to the subfeature bit						
					add_mRNA_to_dict(bits, chr)	
				elif type in ["exon", "CDS"]:
					add_exon_or_cds_to_dict(bits, chr)



def parse_info(info_bit, flag):
	#print(flag, info_bit)
	if flag in "Name":
		regex = "Name=(.*?);"
	elif flag =="ID":	
		regex = "ID=(.*?);"
	elif flag == "parent":
		regex = "Parent=(.*?);"	
	elif flag == "GeneID":
		regex = "=GeneID:(.*?)[,;]"
	elif flag == "Gene":
		regex = "gene=(.*?);"		
	elif flag == "Alias":
		regex = "Alias=(.*);*"	
	re_out = re.search(regex, info_bit)
	
	if re_out is None:
		re_out = [0,"NA"]
	return(re_out[1])

def revcom(seq):
	seq      = seq.upper()
	com_dict = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
	rev      = seq[::-1]
	revcom   = "".join([com_dict[i] for i in rev])
	return(revcom)



#import fasta into seqIO record dictionary:
record_dict = {}
parse_fasta(args.fasta)#record_dict = SeqIO.index(args.fasta, "fasta") 					 #this way is apparently faster for really big files like genomes. 
#print(record_dict)

#update record dictionary wiuth feature stuff"
# the gff dict is at the top level each chr in the fasta - gff_dict keys should be a subset of the fasta dict keys.
#inside is each feature in it's own dict:
#features have a type, source, start, end, something1, strand, frame, info, and subfeature entries. 	
gff_dict = {k:{"region":{}, "gene":{}, "pseudogene":{}, "lnc_RNA":{}, "tRNA":{}, "rRNA":{}, "pseudogenic_transcript":{}, "ncRNA_gene":{}} for k in record_dict} #make the gff_dict have the same keys as the record_dict 
parent_tree_dict = {}
parent_type_dict = {}
alias_dict = {}
parse_gff(args.GFF)
#pprint.pprint(gff_dict["seq2"]["gene"])

#now extract the proper fasta bits:
extract_cds(args.Output_type) #output type can be "fasta" or "GC_table"






	