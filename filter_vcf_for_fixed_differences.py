#!/usr/bin/env python3

import re
import sys

if len(sys.argv)<2:
	expl="\nfilter_vcf_for_fixed_differences.py\ncreated by TDB 11-04-2014\n\nThis script takes one input at the command line:\n\t(1) a vcf file of 14 campbelli and 11 sungorus founders that have been run through the gatk snp calling pipeline\n\nit generates through stdout a vcf file that has only fixed differences in it.\n\n"
	print(expl, file=sys.stderr)
else:
	vcf=open(sys.argv[1],"r")
	for line in vcf:
		line=line.strip("\n")
		if line.startswith("#"):
			print(line, file=sys.stdout)
		else:
			fields=line.split("\t")
			#print(fields[9],fields[28], file=sys.stderr)
			B1=fields[9].split(":")[0]
			B2=fields[10].split(":")[0]
			B3=fields[11].split(":")[0]
			B4=fields[12].split(":")[0]
			B5=fields[13].split(":")[0]
			B6=fields[14].split(":")[0]
			B7=fields[15].split(":")[0]
			B8=fields[16].split(":")[0]
			B9=fields[17].split(":")[0]
			B10=fields[18].split(":")[0]
			S1=fields[19].split(":")[0]
			S2=fields[20].split(":")[0]
			S3=fields[21].split(":")[0]
			S4=fields[22].split(":")[0]
			S5=fields[23].split(":")[0]
			S6=fields[24].split(":")[0]
			S7=fields[25].split(":")[0]
			S8=fields[26].split(":")[0]
			S9=fields[27].split(":")[0]
			S10=fields[28].split(":")[0]
			if(B1==B2==B3==B4==B5==B6==B7==B8==B9==B10):
				if(S1==S2==S3==S4==S5==S6==S7==S8==S9==S10!=B1):
					print(line, file=sys.stdout)




