#Converts a FASTA file to a tab-delimited file, assuming the fasta file is in the working directory
#Change INPUT.FASTA to FASTA file to be converted
# Change OUTPUT.TAB to whatever you want the output name to be

perl -e '
$count=0; 
$len=0; 
while(<>) { 
	s/\r?\n//; 
	s/\t/ /g; 
	if (s/^>//) { 
		if ($. != 1) { 
			print "\n" 
		} 
		s/ |$/\t/; 
		$count++; 
		$_ .= "\t"; 
	} else { 
		s/ //g; 
		$len += length($_) 
	} 
	print $_; 
}
print "\n";
' $1 > $2
