#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Const::Fast;
use Path::Class;
use Smart::Comments;

## Skip codon if reference is het
const my %AMINO_ACIDS => ( 'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I', 'ATG' => 'M',
			   'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
			   'AAC' => 'N', 'AAT' => 'N', 'AAA' => 'K', 'AAG' => 'K',
			   'AGC' => 'S', 'AGT' => 'S', 'AGA' => 'R', 'AGG' => 'R',
			   'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',
			   'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
			   'CAC' => 'H', 'CAT' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
			   'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',
			   'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
			   'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
			   'GAC' => 'D', 'GAT' => 'D', 'GAA' => 'E', 'GAG' => 'E',
			   'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
			   'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',
			   'TTC' => 'F', 'TTT' => 'F', 'TTA' => 'L', 'TTG' => 'L',
			   'TAC' => 'Y', 'TAT' => 'Y', 'TAA' => '*', 'TAG' => '*',
			   'TGC' => 'C', 'TGT' => 'C', 'TGA' => '*', 'TGG' => 'W' );
const my @BASES => qw( A C G T );
const my %EXPAND_IUPAC => ( 'M' => [ 'A', 'C' ], 'R' => [ 'A', 'G' ], 'W' => [ 'A', 'T' ], 'S' => [ 'C', 'G' ], 'Y' => [ 'C', 'T' ], 'K' => [ 'G', 'T' ],
			    'V' => [ 'A', 'C', 'G' ], 'H' => [ 'A', 'C', 'T' ], 'D' => [ 'A', 'G', 'T' ], 'B' => [ 'C', 'G', 'T' ],
			    'N' => [ 'A', 'C', 'G', 'T' ] );
#const my $A_THALIANA_GENOME_FILE => '/Users/azs607/Bangor/afzelia/a_thaliana_ref/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa';
#const my $A_THALIANA_GFF_FILE => '/Users/azs607/Bangor/afzelia/a_thaliana_ref/Arabidopsis_thaliana.TAIR10.37.gff3';

GetOptions(
    "genome|g=s" => \my $A_THALIANA_GENOME_FILE,
    "gff3|a=s" => \my $A_THALIANA_GFF_FILE,
    "vcf|v=s" => \my $vcf_file, # VCF file from alignments of A. lyrata sequencing reads against A. thaliana reference
    "out|o=s" => \my $out_file
    );
die "USAGE:  perl calc_dN_dS.pl -v <vcf_file> -o <out_file> -g <genome file> -a <annotation file (gff3)>\n" unless defined $vcf_file and defined $out_file and defined $A_THALIANA_GENOME_FILE and defined $A_THALIANA_GFF_FILE ;

my $a_thaliana_cdna_coords = get_cdna_coordinates();

my $a_lyrata_variations = parse_vcf_file( $vcf_file );

my $out_fh = file( $out_file )->openw;
$out_fh->print( '#Transcript,No. silent sites,No. replacement site,No. silent changes,No. replacement changes,dN,dS,dNdS' );
for my $transcript( sort keys %$a_thaliana_cdna_coords ){
    ### Transcript: $transcript
    my ( $a_thaliana_cdna_seq, $a_lyrata_cdna_seq ) = get_cdna_sequences( $a_thaliana_cdna_coords->{ $transcript }, $a_lyrata_variations );
    next unless length $a_thaliana_cdna_seq == length $a_lyrata_cdna_seq and length $a_thaliana_cdna_seq > 10; # Can't currently deal with indels
    unless ( $a_thaliana_cdna_seq =~ /^[ACGT]$/ ){ # Can't currently deal with heterozygous bases in reference unless they give same aa
	$a_thaliana_cdna_seq = replace_heterozygous_bases_if_synonymous( $a_thaliana_cdna_seq );
	next unless $a_thaliana_cdna_seq =~ /^[ACTG]+$/; 
	$a_lyrata_cdna_seq = replace_heterozygous_bases_if_synonymous( $a_lyrata_cdna_seq );
	next unless $a_lyrata_cdna_seq =~ /^[ACGT]+$/;
    }

    my ( $orf1 ) = get_orf( $a_thaliana_cdna_seq, $transcript . ' - A. thaliana cDNA' );
    my ( $orf2 ) = get_orf( $a_lyrata_cdna_seq, $transcript . ' - A. lyrata cDNA' );
    my ( $nr_S_in_ORF, $nr_N_in_ORF ) = calc_SN_per_ORF( $orf1 );
    my ( $nr_S_between_seqs, $nr_N_between_seqs ) = calc_SN_between_seqs( $orf1, $orf2, $transcript . ' - A. thaliana cDNA',  $transcript . ' - A. lyrata cDNA' );
    my $dN = $nr_N_between_seqs / $nr_N_in_ORF;
    my $dS = $nr_S_between_seqs / $nr_S_in_ORF;
    my $dNdS = $dS == 0 ? 'undefined' : $dN / $dS;

    $out_fh->print( "\n" . join( ',', $transcript, $nr_S_in_ORF, $nr_N_in_ORF, $nr_S_between_seqs, $nr_N_between_seqs, $dN, $dS, $dNdS ) );
}

sub calc_SN_between_seqs{
    my ( $orf1, $orf2, $id1, $id2 ) = @_;
    
    die "ORFs for $id1 and $id2 are unequal in length\n" unless length $orf1 == length $orf2;
    my $total_S_between_seqs = 0;
    my $total_N_between_seqs = 0;
    my $codon_start_ix = 0;
    while( $codon_start_ix  < length $orf1 ){
	my ( $codon_S_between_seqs, $codon_N_between_seqs ) = compare_codons( substr( $orf1, $codon_start_ix, 3 ), substr( $orf2, $codon_start_ix, 3 ) );
	$total_S_between_seqs += $codon_S_between_seqs;
	$total_N_between_seqs += $codon_N_between_seqs;
	$codon_start_ix += 3;
    }

    return ( $total_S_between_seqs, $total_N_between_seqs );
}
	
sub compare_codons{
    my ( $codon1, $codon2 ) = @_;

    my $nr_differences = ( $codon1 ^ $codon2 ) =~ tr/\0//c;

    return ( 0, 0 ) if $nr_differences == 0;
    return get_SN_per_site( $codon1, $codon2 ) if $nr_differences == 1;

    if( $nr_differences == 2 ){
	my ( $int_codon1, $int_codon2 ) = get_intermediate_codons_for_2_differences( $codon1, $codon2 );
	my ( $route1_step1_S, $route1_step1_N ) = get_SN_per_site( $codon1, $int_codon1 );
	my ( $route1_step2_S, $route1_step2_N ) = get_SN_per_site( $int_codon1, $codon2 );
	my ( $route2_step1_S, $route2_step1_N ) = get_SN_per_site( $codon1, $int_codon2 );
	my ( $route2_step2_S, $route2_step2_N ) = get_SN_per_site( $int_codon2, $codon2 );
	return ( mean( $route1_step1_S, $route2_step2_S ) + mean( $route1_step2_S, $route2_step1_S ),
		 mean( $route1_step1_N, $route2_step2_N ) + mean(  $route1_step2_N, $route2_step1_N ) );
    }

    die "Something has gone horribly wrong - >3 differences between codons $codon1 and $codon2 !!!\n" unless $nr_differences == 3;
    my ( $int1_codon1, $int1_codon2, $int1_codon3 ) = get_first_step_intermediate_codons_for_3_differences( $codon1, $codon2 );
    my ( $int2_codon1, $int2_codon2 ) = get_intermediate_codons_for_2_differences( $int1_codon1, $codon2 );
    my ( $int2_codon3, $int2_codon4 ) = get_intermediate_codons_for_2_differences( $int1_codon2, $codon2 );
    my ( $int2_codon5, $int2_codon6 ) = get_intermediate_codons_for_2_differences( $int1_codon3, $codon2 );
    my ( $route12_step1_S, $route12_step1_N ) = get_SN_per_site( $codon1, $int1_codon1 );
    my ( $route34_step1_S, $route34_step1_N ) = get_SN_per_site( $codon1, $int1_codon2 );
    my ( $route56_step1_S, $route56_step1_N ) = get_SN_per_site( $codon1, $int1_codon3 );
    my ( $route1_step2_S, $route1_step2_N ) = get_SN_per_site( $int1_codon1, $int2_codon1 );
    my ( $route2_step2_S, $route2_step2_N ) = get_SN_per_site( $int1_codon1, $int2_codon2 );
    my ( $route3_step2_S, $route3_step2_N ) = get_SN_per_site( $int1_codon2, $int2_codon3 );
    my ( $route4_step2_S, $route4_step2_N ) = get_SN_per_site( $int1_codon2, $int2_codon4 );
    my ( $route5_step2_S, $route5_step2_N ) = get_SN_per_site( $int1_codon3, $int2_codon5 );
    my ( $route6_step2_S, $route6_step2_N ) = get_SN_per_site( $int1_codon3, $int2_codon6 );
    my ( $route1_step3_S, $route1_step3_N ) = get_SN_per_site( $int2_codon1, $codon2 );
    my ( $route2_step3_S, $route2_step3_N ) = get_SN_per_site( $int2_codon2, $codon2 );
    my ( $route3_step3_S, $route3_step3_N ) = get_SN_per_site( $int2_codon3, $codon2 );
    my ( $route4_step3_S, $route4_step3_N ) = get_SN_per_site( $int2_codon4, $codon2 );
    my ( $route5_step3_S, $route5_step3_N ) = get_SN_per_site( $int2_codon5, $codon2 );
    my ( $route6_step3_S, $route6_step3_N ) = get_SN_per_site( $int2_codon6, $codon2 );
    return( mean( $route12_step1_S, $route34_step1_S, $route56_step1_S )
	    + mean( $route1_step2_S, $route2_step2_S, $route3_step2_S, $route4_step2_S, $route5_step2_S, $route6_step2_S )
	    + mean( $route1_step3_S, $route2_step3_S, $route3_step3_S, $route4_step3_S, $route5_step3_S, $route6_step3_S ),
	    mean( $route12_step1_S, $route34_step1_S, $route56_step1_S )
	    + mean( $route1_step2_N, $route2_step2_N, $route3_step2_N, $route4_step2_N, $route5_step2_N, $route6_step2_N )
	    + mean( $route1_step3_N, $route2_step3_N, $route3_step3_N, $route4_step3_N, $route5_step3_N, $route6_step3_N ) );
}

sub calc_SN_per_codon{
    my $codon = shift;

    my $nr_S = 0;
    my $nr_N = 0;
    
    for my $codon_ix( 0 .. 2 ){
	for my $base( @BASES ){
	    my $mutated_base = substr( $codon, $codon_ix, 1 );
	    next if $mutated_base eq $base;
	    my $before_mutated_base = $codon_ix == 0 ? '' : substr( $codon, 0, $codon_ix );
	    my $after_mutated_base = $codon_ix == 2 ? '' : substr( $codon, $codon_ix + 1 );
	    my $mutated_codon = $before_mutated_base . $base . $after_mutated_base;
	    if ( $AMINO_ACIDS{ $codon } eq $AMINO_ACIDS{ $mutated_codon } ){
		$nr_S++;
	    }
	    else{
		$nr_N++;
	    }
	}
    }
    die "calc_SN_per_codon error: S + N != 3\nS = $nr_S, N=$nr_N\n" unless ( $nr_S + $nr_N ) / 3 == 3;

    return ( $nr_S / 3, $nr_N / 3 );
}

sub calc_SN_per_ORF{
    my $orf = shift;

    my $total_S = 0;
    my $total_N = 0;

    my $codon_start_ix = 0;
    while( $codon_start_ix  < length $orf ){
	my $codon = substr( $orf, $codon_start_ix, 3 );
	my ( $codon_S, $codon_N ) = calc_SN_per_codon( $codon );
	$total_S += $codon_S;
	$total_N += $codon_N;
	$codon_start_ix += 3;
    }

    return ( $total_S, $total_N );
}

sub get_cdna_coordinates{
    my %cdna_coord;
    my $fh = file( $A_THALIANA_GFF_FILE )->openr;
    while( my $line = $fh->getline() ){
	next if $line =~ /^#/;
	chomp $line;
	my ( $chr, $source, $feature, $start, $end, $score, $strand, $frame, $attributes ) = split( "\t", $line );
	next unless $feature eq 'CDS';
	my ( $parent_transcript ) = $attributes =~ /Parent=transcript:(AT.G\d+\.\d+);/;
	$cdna_coord{ $parent_transcript }{ 'chromosome' } = $chr;
	$cdna_coord{ $parent_transcript }{ 'strand' } = $strand;
	$cdna_coord{ $parent_transcript }{ 'coordinates' }{ $start } = $end;
    }

    return \%cdna_coord;
}

sub get_cdna_sequences{
    my ( $transcript_coords, $variations ) = @_;
    
    my $a_thaliana_sequence = '';
    my $a_lyrata_sequence = '';
    for my $start( sort {$a<=>$b} keys $transcript_coords->{ 'coordinates' } ){
	my $chr = $transcript_coords->{ 'chromosome' };
	my $end = $transcript_coords->{ 'coordinates' }{ $start };
 
	my $variations_within_cds = get_relative_variation_positions( $start, $end, $variations->{$chr} );
	my $a_thaliana_cds_seq = '';
	open( FAIDX, "samtools faidx $A_THALIANA_GENOME_FILE $chr:$start-$end |" );
	while( <FAIDX> ){
	    next if $_ =~ /^>/;
	    chomp;
	    $a_thaliana_cds_seq .= uc( $_ );
	}
	close( FAIDX );
	$a_thaliana_sequence .= $a_thaliana_cds_seq;
	$a_lyrata_sequence .= incorporate_variations( $a_thaliana_cds_seq, $variations_within_cds );
    }

    return ( $a_thaliana_sequence, $a_lyrata_sequence ) if $transcript_coords->{ 'strand' } eq '+';
    return ( reverse_complement( $a_thaliana_sequence ), reverse_complement( $a_lyrata_sequence ) );
}
		   
sub get_first_step_intermediate_codons_for_3_differences{
    my ( $codon1, $codon2 ) = @_;

    my @codon1_bases = split( //, $codon1 );
    my @codon2_bases = split( //, $codon2 );

    my $int_codon1 = $codon2_bases[0] . $codon1_bases[1] . $codon1_bases[2];
    my $int_codon2 = $codon1_bases[0] . $codon2_bases[1] . $codon1_bases[2];
    my $int_codon3 = $codon1_bases[0] . $codon1_bases[1] . $codon2_bases[2];

    return( $int_codon1, $int_codon2, $int_codon3 );
}

sub get_intermediate_codons_for_2_differences{
    my ( $codon1, $codon2 ) = @_;
 
    my @codon1_bases = split( //, $codon1 );
    my @codon2_bases = split( //, $codon2 );

    my ( $int_codon1, $int_codon2 );
    if ( $codon1_bases[0] eq $codon2_bases[0] ){
	$int_codon1 = $codon1_bases[0] . $codon1_bases[1] . $codon2_bases[2];
	$int_codon2 = $codon1_bases[0] . $codon2_bases[1] . $codon1_bases[2];
    }
    elsif( $codon1_bases[1] eq $codon2_bases[1] ){
	$int_codon1 = $codon1_bases[0] . $codon1_bases[1] . $codon2_bases[2];
	$int_codon2 = $codon2_bases[0] . $codon1_bases[1] . $codon1_bases[2];
    }
    else{
	die "No matching bases between codons $codon1 and $codon2 when 1 match expected\n" unless $codon1_bases[2] eq $codon2_bases[2];
	$int_codon1 = $codon2_bases[0] . $codon1_bases[1] . $codon1_bases[2];
	$int_codon2 = $codon1_bases[0] . $codon2_bases[1] . $codon1_bases[2];
    }

    return ( $int_codon1, $int_codon2 );
}

sub get_orf{
    my ( $orf, $id ) = @_;

    my $flag = 0;
    $orf = uc( $orf );
    unless ( $orf =~ /^ATG/ ){
	$flag++;
	warn( "No start codon for $id\n" );
    }
    unless ( $orf =~ /TA[AG]$/ or $orf =~ /TGA$/ ){
	$flag += 2;
	warn( "No stop codon for $id\n" );
    }
    my $part_codon_bases = ( length $orf ) % 3;
    if( $part_codon_bases != 0 ){
	$flag += 4;
	warn( "Length of $id ORF is not a multiple of 3\n" );
	if ( $flag == 5 ){
	    warn "Trimming bases from start of $id ORF to keep in frame\n";
	    return substr( $orf, $part_codon_bases );
	}
	else{
	    if ( $flag == 6 ){
		warn "Trimming bases from end of $id ORF to keep in frame\n";
	    }
	    else{
		warn( "Cannot determine correct way to trim $id ORF - trimming bases from end to keep in frame\n" );
	    }
	    return substr( $orf, 0, ( length $orf ) - $part_codon_bases );
	}
    }

    return $orf;
}

sub get_relative_variation_positions{
    my ( $start, $end, $chromosome_variations ) = @_;
    my %relative_positions;
    return \%relative_positions unless defined $chromosome_variations;

    for my $pos( sort {$a<=>$b} %$chromosome_variations ){
	last if $pos > $end;
	next if $pos < $start;
	my $relative_pos = $pos - $start;
	$relative_positions{ $relative_pos }{ 'ref' } = $chromosome_variations->{ $pos }{ 'ref' };
	$relative_positions{ $relative_pos }{ 'alt' } = $chromosome_variations->{ $pos }{ 'alt' };
    }

    return \%relative_positions;
}

sub get_SN_per_site{
    my ( $codon1, $codon2 ) = @_;

    return ( 1, 0 ) if $AMINO_ACIDS{ $codon1 } eq $AMINO_ACIDS{ $codon2 };
    return ( 0, 1 );
}

sub incorporate_variations{
    my ( $seq, $variations ) = @_;

    my $indel_offset = 0;
    for my $pos( sort {$a<=>$b} keys %$variations ){
	my $seq_before_var = substr( $seq, 0, $pos + $indel_offset );
	my $seq_to_change = substr( $seq, $pos + $indel_offset, length( $variations->{ $pos }{ 'ref' } ) );
	my $seq_after_var = substr( $seq, $pos + $indel_offset + length( $variations->{ $pos }{ 'ref' } ) );
	$seq = $seq_before_var . $variations->{ $pos }{ 'alt' } . $seq_after_var;
	$indel_offset += length( $variations->{ $pos }{ 'alt' } ) - length( $variations->{ $pos }{ 'ref' } );
    }

    return $seq;
}

sub mean{
    my @array = @_;

    my $sum = 0;
    for my $element( @array ){
	$sum += $element;
    }

    return $sum / ( scalar @array );
}

sub parse_vcf_file{
    my $file = shift;

    my %variations;
    my $fh = file( $file )->openr;
    while( my $line = $fh->getline() ){
	next if $line =~ /^#/;
	chomp $line;
	my ( $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $attributes) = split( "\t", $line );
	next if $alt =~ /,/;
	#if ( defined $attributes ){
	#    my ( $genotype ) = $attributes =~ /^(\d\/\d):/;
	#    next unless $genotype eq '1/1'; #only interested in homozygous variations
	#}
	$variations{ $chr }{ $pos }{ 'ref' } = $ref;
	$variations{ $chr }{ $pos }{ 'alt' } = $alt;
    }

    return \%variations;
}

sub replace_heterozygous_bases_if_synonymous{
    my $sequence = shift;

    for my $base_ix( 0 .. ( length $sequence ) - 1 ){
	my $base = substr( $sequence, $base_ix, 1 );
	next if $base =~ /^[ACGT]$/;
	my @possible_bases = $EXPAND_IUPAC{ $base };
	my $codon_modulus = $base_ix % 3;
	my $codon_before_iupac = substr( $sequence, $base_ix - $codon_modulus, $codon_modulus );
	my $codon_after_iupac = substr( $sequence, $base_ix + 1, 2 - $codon_modulus );
	my @possible_codons;
	for my $possible_base( @possible_bases ){
	    push @possible_codons, $codon_before_iupac . $possible_base . $codon_after_iupac;
	}
	for my $possible_codon( @possible_codons ){
	    return $sequence unless $possible_codon =~ /^[ACGT]{3}$/;
	    return $sequence unless $AMINO_ACIDS{ $possible_codon } eq $AMINO_ACIDS{ $possible_codons[0] };
	}
	$sequence = substr( $sequence, 0, $base_ix ) . $possible_bases[0] . substr( $sequence, $base_ix + 1 );
    }

    return $sequence;
}

sub reverse_complement{
    my $seq = shift;

    $seq =~ s/A/t/g;
    $seq =~ s/T/a/g;
    $seq =~ s/C/g/g;
    $seq =~ s/G/c/g;

    return uc( reverse( $seq ) );
}
