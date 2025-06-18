#!/usr/bin/perl -w
############################################################################
# Module:      Sequence 
# Description: Sequence analysis functions
# History:     Rob Gifford, Decemeber 2006: Creation
############################################################################
package Sequence;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

############################################################################
# Globals
############################################################################
1;

# Create hash with IUPAC codes for nucleic acids
my %iupac;
$iupac{A} = 'adenine';
$iupac{T} = 'thymine';
$iupac{G} = 'guanine';
$iupac{C} = 'cytosine';

# Twofold degenerate nucleotides
my @R = qw [ G A ];
$iupac{R} = \@R;
my @Y = qw [ T C ];
$iupac{Y} = \@Y;
my @K = qw [ G T ];
$iupac{K} = \@K;
my @M = qw [ A C ];
$iupac{M} = \@M;
my @S = qw [ G C ];
$iupac{S} = \@S;
my @W = qw [ A T ];
$iupac{W} = \@W;

# Threefold degenerate nucleotides
my @B = qw [ G T C ];
$iupac{B} = \@B;
my @D = qw [ G A T ];
$iupac{D} = \@D;
my @H = qw [ A C T ];
$iupac{H} = \@H;
my @V = qw [ G C A ];
$iupac{V} = 1;

# Fourfold degenerate nucleotides
my @N = qw [ A G C T ];
$iupac{N} = \@N;

my @dinucs = qw [ AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT ];

# Create base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create a new Sequence.pm object
#***************************************************************************
sub new {

	my ($invocant, $sequence, $header, $id) = @_;
	
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
		iupac         => \%iupac,
		dinucs        => \@dinucs,
		sequence      => $sequence,
		sequence_id   => $id,
		header        => $header,
	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Manipulations of sequence string
############################################################################

#***************************************************************************
# Subroutine:  translate
# Description: Translate na sequence to aa
# Arguments:   $sequence: sequence to translate
# Arguments:   $frame to translate in: (possible values 1, 2, 3)
#***************************************************************************
sub translate {

	my ($self, $sequence, $frame) = @_;

	
	my @sequence = split ('', $sequence);
	
	# do frame adjustment
	#unless ($frame) { die; }
	unless ($frame) { $frame = 1; }
	if    ($frame eq 1)  {  }
	elsif ($frame eq 2)  { shift (@sequence); }
	elsif ($frame eq 3)  { shift (@sequence);	
		                   shift (@sequence); }
	else {
		die "\n\t translation in frame '$frame' not implemented\n\n";
	}

	# Translate the sequence
	my $codon;
	my $translation;
	my $codon_index;
	foreach my $nt (@sequence) {
		$codon .= $nt;
		$codon_index++;
		if ($codon_index eq 3) {
			if ($codon =~ /-/) {
				$translation .= '-';
			}
			# can't do degenerate codons
			elsif ($codon =~ 'N') {
				$translation .= '?';
			}
			else {
				my @codon_list;
				$self->get_codon_list($codon, \@codon_list);
				# if the codon is degenerate arbitrarily pick first one in the list
				my $list_codon = $codon_list[0];
				my $aa = $self->codon2aa($list_codon);
				$translation .= $aa;
			}
			$codon = '';
			$codon_index = 0;
		}
	}
	return $translation;
}

#***************************************************************************
# Subroutine:  reverse and complement
# Description: reverse and complement sequence
# Arguments:   $sequence: the sequence string to reverese and complement 
#***************************************************************************
sub reverse_and_complement {

	my ($self, $sequence) = @_;
    $sequence = uc $sequence;
	
	my @sequence = split ('', $sequence);

	@sequence = reverse(@sequence);

	my $rc_sequence;
	foreach my $nt (@sequence) {
		
		if    ($nt eq 'A') { $nt = 'T'; }
		elsif ($nt eq 'T') { $nt = 'A'; }
		elsif ($nt eq 'G') { $nt = 'C'; }
		elsif ($nt eq 'C') { $nt = 'G'; }
		elsif ($nt eq 'N') { $nt = 'N'  }
	
		else {
			print "\n\t ERROR: input was $sequence";
			die "\n\t This routine cannot yet deal with NT $nt\n"; 
		}
	
		$rc_sequence .= $nt;
	}
	return $rc_sequence;
}

#***************************************************************************
# Subroutine:  extract_subsequence
# Description: extract subseqeunce region from a given sequence 
# Arguments:   $sequence: the sequence string to extract from 
#              $start: region start position
#              $stop:  region start position
#***************************************************************************
sub extract_subsequence {
	
	my ($self, $sequence, $start, $stop) = @_;

	# Sanity checking
	unless ($start) { $start = 1; }
	unless ($sequence and $start and $stop) { 
		die;
	}
	
	# Do the extraction
	my @chars = split ('', $sequence);
	my $subsequence = '';
	my $i = 0;
	foreach my $char (@chars) {
		$i++;
		if ($i >= $start and $i <= $stop) {
			$subsequence .= $char;
		}
	}
	return $subsequence;
}

#***************************************************************************
# Subroutine:  get_translations 
# Description: Get all possible translations of a given codon
# Arguments:   $data_ref: reference to hash to store all distinct AAs that 
#              can be obtained by translating the (possibly degenerate) 
#              codon
#***************************************************************************
sub get_translations {

	my ($self, $codon) = @_;

	# translate the codon (to all possible if degenerate)
	my @codon_list;
	$self->get_codon_list($codon, \@codon_list);
	my $translations = undef;
	my %seen_translations;
	foreach my $possible_codon (@codon_list) {
	
		my $translation .= $self->codon2aa($possible_codon);
		unless ($seen_translations{$translation}) {
			$translations .= $translation;
			$seen_translations{$translation} = 1;
		}
	}

	return $translations;
}

#***************************************************************************
# Subroutine:  get_codon_list
# Description: get a list with all the possible codons from a degenerate 
#              codon
# Arguments:   $codon: the degenerate codon to deal with
#              $array_ref: reference to an array to store the codon list
#***************************************************************************
sub get_codon_list {

	my($self, $codon, $array_ref) = @_;

	# populate a hash of arrays with degeneracy codes
    my (%degeneracy_code) = (
    
		'A' => ['A'],
		'T' => ['T'],
		'C' => ['C'],
		'G' => ['G'],
		'-' => ['-'],
		'M' => ['A' , 'C'],
		'R' => ['A' , 'G'],
		'W' => ['A' , 'T'],
		'S' => ['G' , 'C'],
		'Y' => ['C' , 'T'],
		'K' => ['G' , 'T'],
		'H' => ['A' , 'C', 'T'],
		'V' => ['A' , 'G', 'C'],
		'D' => ['A' , 'G', 'T'],
		'B' => ['G' , 'T', 'C'],
		'N' => ['A' , 'G', 'C' , 'T'],
	);
	
	my @bases = split('', $codon);
	my $first  = $degeneracy_code{$bases[0]};
	my $second = $degeneracy_code{$bases[1]};
	my $third  = $degeneracy_code{$bases[2]};

	foreach my $first_pos (@$first) {

		foreach my $second_pos (@$second) {
		
			foreach my $third_pos (@$third) {
		
				my $codon_option = $first_pos . $second_pos . $third_pos;
				push (@$array_ref, $codon_option);
			}
		}
	}
}

#***************************************************************************
# Subroutine:  codon2aa
# Description: does translation accorsing to standard genetic code
#              (taken from the o'reilly book 'beginning perl for 
#              bioinformatics).
#***************************************************************************
sub codon2aa {
   
	my ($self, $codon) = @_;
    
	$codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }
	else {
		if ($codon =~ /[-~MRWSYKHVDBN]/) {  
			#print "\n\t Can't translate codon $codon";
			return '-';
		}	
		else {
			print "Bad codon '$codon'!!\n";
			return '-';
		}
	}
}

#***************************************************************************
# Subroutine:  index_sequence 
# Description: 
#***************************************************************************
sub index_sequence {

	my ($self, $sequence, $indexed_sequence_ref, $start) = @_;
	
	my @sequence = split ('', $sequence);
	my $index = 0;
	if ($start) { $index = $start; }
	foreach my $residue (@sequence) {
		$index++;
		#print "\n\t INDEX $index $residue";
		$indexed_sequence_ref->{$index} = $residue;
	}
}

############################################################################
# Codon level operations (fxns that act on codons)
############################################################################

#***************************************************************************
# Subroutine:  compare_codons 
# Description: 
#***************************************************************************
sub compare_codons {

	my ($self, $ref_codon, $codon, $codon_data) = @_;
	
	# Matrix
	my @matrix;
	push (@matrix, 'A->T');
	push (@matrix, 'A->C');
	push (@matrix, 'A->G');
	push (@matrix, 'T->A');
	push (@matrix, 'T->C');
	push (@matrix, 'T->G');
	push (@matrix, 'C->A');
	push (@matrix, 'C->T');
	push (@matrix, 'C->G');
	push (@matrix, 'G->A');
	push (@matrix, 'G->T');
	push (@matrix, 'G->C');
	foreach my $change (@matrix) {
		$codon_data->{$change} = 0;
	}

	# Create variables 
	my $num_syn = 0;
	my $num_nonsyn = 0;
	my $num_changes = 0;
	my @codon     = split('', $codon);
	my @ref_codon = split('', $ref_codon);
	my $ref_aa    = $self->codon2aa($ref_codon);
	my $i = 0;
	$codon_data->{num_transitions} = 0;
	$codon_data->{num_transversions} = 0;
	foreach my $nt (@codon) {

		# get the reference nt
		my $ref_nt = $ref_codon[$i];
		unless ($ref_nt eq $nt) { 
			
			$num_changes++;
			
			# Deal with any degeneracy
			my @site;
			$self->decompose_base_mixture($nt, \@site);
			my $degeneracy = scalar @site;
			unless ($degeneracy > 0) { die; }
			my $increment = 1 / $degeneracy;
			
			# Carry out nuc vs nuc tests
			foreach my $d_nt (@site) {	
				
				# Work out transitions/transversions
				#print "\n\t ref_nt $ref_nt and $d_nt and $increment and $codon_data\n\n";
				unless ($ref_nt and $d_nt and $increment and $codon_data) { die; }
				$self->compare_nucs($ref_nt, $d_nt, $increment, $codon_data);
				
				# Work out if this is a syn or non-syn change
				my @test_codon = @ref_codon;
				$test_codon[$i] = $d_nt;
				my $test_codon = join ('', @test_codon);
				my $aa = $self->codon2aa($test_codon); 
				if ($aa eq $ref_aa) {
					#$devtools->print_array(\@site); 
					#$devtools->print_hash($codon_data); 
					#die;
					$num_syn = $num_syn + $increment;	
				}
				else {
					$num_nonsyn = $num_nonsyn + $increment;	
				}
			
				# Populate change matrix
				my $change_key = $ref_nt . '->' . $d_nt;	
				if ($codon_data->{$change_key}) { 
					$codon_data->{$change_key} = $codon_data->{$change_key} + $increment;
				}
				else {
					$codon_data->{$change_key} = $increment;
				}
			}
		}

		# Increment within codon site count
		$i++;
	}
	
	$codon_data->{num_changes}  = $num_changes;
	$codon_data->{num_nonsyn}   = $num_nonsyn; 
	$codon_data->{num_syn}      = $num_syn;
	
	#print "\n\t REF $ref_codon vs $codon";
	#$devtools->print_hash($codon_data);
}

#***************************************************************************
# Subroutine:  validate_codon
# Description: test wether a three-char triplet is a valid codon
#***************************************************************************
sub validate_codon {

	my ($self, $codon) = @_;
	my $valid = 1;
	if  ($codon =~ /-/ or $codon =~ / / or $codon eq 'NNN') {  
		$valid = undef;
	}
	return $valid;	
}

#***************************************************************************
# Subroutine:  decompose_base_mixture
# Description: convert degenerate base to a list of possibilities
#***************************************************************************
sub decompose_base_mixture {

	my ($self, $base, $nt_ref) = @_;

	# Not degenerate
	if    ($base eq 'A') { push(@$nt_ref, $base); }
	elsif ($base eq 'T') { push(@$nt_ref, $base); }
	elsif ($base eq 'C') { push(@$nt_ref, $base); }
	elsif ($base eq 'G') { push(@$nt_ref, $base); }

	# Two-fold degenerate
	elsif ($base eq 'M') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'C');
	}
	elsif ($base eq 'R') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'G');
	}
	elsif ($base eq 'W') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'T');
	}
	elsif ($base eq 'S') { 
		push(@$nt_ref, 'G');
		push(@$nt_ref, 'C');
	}
	elsif ($base eq 'Y') { 
		push(@$nt_ref, 'C');
		push(@$nt_ref, 'T');
	}
	elsif ($base eq 'K') { 
		push(@$nt_ref, 'G');
		push(@$nt_ref, 'T');
	}
	
	# Three-fold degenerate
	elsif ($base eq 'H') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'C');
		push(@$nt_ref, 'T');
	}
	elsif ($base eq 'V') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'G');
		push(@$nt_ref, 'C');
	}
	elsif ($base eq 'D') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'G');
		push(@$nt_ref, 'T');
	}
	elsif ($base eq 'B') { 
		push(@$nt_ref, 'G');
		push(@$nt_ref, 'T');
		push(@$nt_ref, 'C');
	}
	
	# Four-fold degenerate
	elsif ($base eq 'N') { 
		push(@$nt_ref, 'A');
		push(@$nt_ref, 'G');
		push(@$nt_ref, 'C');
		push(@$nt_ref, 'T');
	}
	else {
		push(@$nt_ref, $base);
		#print "\n\t Something wrong here: '$base' is not a recognized symbol\n";
	}
}

#***************************************************************************
# Subroutine:  compare_nucs
# Description: 
#***************************************************************************
sub compare_nucs {

	my ($self, $ref_nt, $nt, $increment, $codon_data) = @_;
	
	# Create variables 
	my $num_transitions = 0;
	my $num_transversions = 0;
	
	# Purines
	my %purines;
	$purines{A} = 1;
	$purines{G} = 1;
	my %pyrimidines;
	$purines{C} = 1;
	$purines{T} = 1;

	if (    $purines{$ref_nt} and $purines{$nt} 
	or  $pyrimidines{$ref_nt} and $pyrimidines{$nt} ) {
		$num_transitions = $num_transitions + $increment;
	}
	elsif ( $pyrimidines{$ref_nt} and $purines{$nt} 
	or  $purines{$ref_nt} and $pyrimidines{$nt} ) {
		$num_transversions++;
		$num_transversions = $num_transversions + $increment;
	}

	$codon_data->{num_transversions} = $codon_data->{num_transversions} + $num_transversions;
	$codon_data->{num_transitions} = $codon_data->{num_transitions} + $num_transitions;

}


#***************************************************************************
# Subroutine:  get_frame 
# Description: get reading frame (deprecated possibly stupid)
#***************************************************************************
sub get_frame {

	my ($self, $start) = @_;

	my $frame = 0;
	my $remainder = $start % 3;
	print "\n\t remainder: $start: '$remainder'";
	
	if ($remainder eq 0) { $frame = 3; }
	if ($remainder eq 1) { $frame = 1; }
	if ($remainder eq 2) { $frame = 2; }
	
	return $frame;
}

#***************************************************************************
# Subroutine:  truncated_header
# Description: shorten sequence *header* - old version is stored
#***************************************************************************
sub truncate_header {

	my ($self, $truncated, $accession_from_gb) = @_;

	my $header   = $self->{header};	
	my $sequence = $self->{sequence};

	my $n_header;
	$accession_from_gb = 1;
	if ($accession_from_gb) {

		#my @header  = split('\.', $header);
		my @header  = split('_', $header);
		$n_header = pop @header;
		#$n_header = shift @header;
		#print "\n\t n_header: $n_header"; #die;
		
		unless ($n_header) {
			my $i = 0;
			do {
				my $char = $header[$i];
				if ($char) { $n_header .= $char; }
				$i++;
			} until ($i eq 15);
		}
	}
	else {
		$n_header = $header;
	}

	print "\n\t $header to $n_header"; #die;

	#my @header  = split('.', $header);
	#$n_header = $header[0];
	$n_header =~ s/^\s+//g;
	#$n_header =~ s/\./-/g;
	$n_header =~ s/,/-/g;
	$n_header =~ s/\//-/g;
	$n_header =~ s/\|/-/g;
	$n_header =~ s/\(/-/g;
	$n_header =~ s/\)/-/g;
	$n_header =~ s/\s+/-/g;
	$n_header =~ s/--/-/g;
	my  $fasta = ">$n_header\n$sequence\n";
	push (@$truncated, $fasta);
	$self->{header}     = $n_header;
	$self->{old_header} = $header;
	#$devtools->print_hash($seq_ref);
	#die;
}

############################################################################
# END OF FILE
############################################################################
