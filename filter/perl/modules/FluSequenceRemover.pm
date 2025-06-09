package FluSequenceRemover;

###############################################################################
# Module:       FluSequenceRemover.pm
# Description:  Separates influenza isolate records into those containing
#               reference sequence accessions and those that do not.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(remove_sequences_from_source_directory);

# ------------------ Subroutine Definitions ------------------

sub remove_sequences_from_source_directory {
    my ($datafile1, $datafile2, $fileio, $devtools) = @_;

	my ($datafile1, $datafile2) = @_;

	# The reference sequence data
	my @datafile1;
	$fileio->read_file($datafile1, \@datafile1);
	my $num_lines1 = scalar @datafile1;
	unless ($num_lines1) {
		print "\n\t # No data read from file '$datafile1'\n\n";
		next;		
	}	
	
	# The strain data
	my @datafile2;
	$fileio->read_file($datafile2, \@datafile2);
	my $num_lines2 = scalar @datafile2;
	unless ($num_lines2) {
		print "\n\t # No data read from file '$datafile2'\n\n";
		next;		
	}

	# Index the reference sequence data
	my %refseqIds;
	my $header_line1 = shift @datafile1;
	chomp $header_line1;
	my @header_line1 = split(/\t/, $header_line1);
	foreach my $line (@datafile1) {
		
		chomp $line;
		#print $line;
		my @line = split(/\t/, $line);
		my $source = $line[0];
		my $seqID  = $line[1];
		#print "\n\t source: $source";
		#print "\n\t seqID: $seqID";
		$refseqIds{$seqID} = $source;
		
	}
	#$devtools->print_hash(\%refseqIds); die;
	
	# Iterate through strain data
	my $header_line2 = shift @datafile2;
	chomp $header_line2;
	my @header_line2 = split(/\t/, $header_line2);
	my $i;
	my %header_fields;
	foreach my $field (@header_line2) {
		$i++;
		$header_fields{$i} = $field;
		
	}

	#$devtools->print_hash(\%header_fields); die;
	my @refseq_isolates;
	my @nonrefseq_isolates;
	
	my $k;
	foreach my $line (@datafile2) {
	
	    $k++;
		chomp $line;
		#print "\n\t LINE $k: $line";
		my @line = split(/\t/, $line);
		#$devtools->print_array(\@line);
		my %data;
		my $j;
		foreach my $value (@line) {
		
			$j++;
			my $field = $header_fields{$j};
			$data{$field} = $value;
		
		}
		
		#$devtools->print_hash(\%data); #die;
		my @segments2 = qw [ 1 2 3 4 5 6 7 ];
		my $is_reference = undef;
		foreach my $segment (@segments2) {
			
			$devtools->print_hash(\%data); #die;

			my $key = 'segment' . $segment . '_accession';
			my $seqID = $data{$key};
			#print "\n\t SEQID $seqID  yeah"; die;
			
			if ($refseqIds{$seqID}) {			
				 $is_reference = 'true'; #die;		
			} 
			
		}

		if ($is_reference) {
			#die $line;
			push (@refseq_isolates, "$line\n");
		}
		else {
			push (@nonrefseq_isolates, "$line\n")
		
		}


	}
	
	#$devtools->print_array(\@refseq_isolates);
	unshift(@refseq_isolates, "$header_line2\n");
	my $refseqs_file_name = $datafile2 . '.refseqs.filtered.tsv';
	$fileio->write_file($refseqs_file_name, \@refseq_isolates);
	
	
	#$devtools->print_array(\@nonrefseq_isolates);

	unshift(@nonrefseq_isolates, "$header_line2\n");
	my $nonrefseqs_file_name = $datafile2 . '.nonrefseqs.filtered.tsv';
	$fileio->write_file($nonrefseqs_file_name, \@nonrefseq_isolates);

}

1;  # End of FluSequenceRemover module
