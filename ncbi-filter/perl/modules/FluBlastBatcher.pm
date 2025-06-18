#!/usr/bin/perl -w
############################################################################
# Module:      FluBlastBatcher.pm 
# Description: 
# History:     Rob Gifford January 2024: Creation
############################################################################
package FluBlastBatcher;

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

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new FluBlastBatcher.pm 'object'
#***************************************************************************
sub new {

    my ($class, $fileio, $seqio, $blast_obj, $devtools) = @_;

	# Member variables
    my $self = {
        fileio    => $fileio,
        seqio     => $seqio,
        blast_obj => $blast_obj,
        devtools  => $devtools,
    };
    bless $self, $class;
    return $self;
}

############################################################################
# Top level fxns
############################################################################

#***************************************************************************
# Subroutine:  batch_blast_against_rsl
# Description: Batch BLAST against reference sequence library to find most
#              similar sequences to those in an input FASTA file
#***************************************************************************
sub batch_blast_against_rsl {

	my ($self, $infile) = @_;

	unless ($lib_path and $tmp_path and $output_path) { die; }

	# Read FASTA input
	my @input_fasta;
    $self->{seqio}->read_fasta($infile, \@input_fasta);
    
	# BLAST each sequence in the input file against the RSL	
	my %not_found;
	my %accessions;
	my %full_record;
	foreach my $seq_ref (@input_fasta) {

		# Get required data about the query sequence
		my $seq_id   = $seq_ref->{sequence_id};
		my $header   = $seq_ref->{header};
		my $sequence = $seq_ref->{sequence};

		# Make a FASTA query file
		my $fasta      = ">$seq_id\n$sequence";
		my $query_file = $tmp_path . '/TEMPORARY.fna';
		$self->{fileio}->write_text_to_file($query_file, $fasta);
		my $result_file = $tmp_path . "/$seq_id.blast_result.tsv";

		# Execute the call to BLAST and parse the results
		$blast_obj->blast('/Users/rob/blast/bin/blastn', $lib_path, $query_file, $result_file);
		my @results;
		$self->{blast_obj}->parse_tab_format_results($result_file, \@results);
		#$devtools->print_array(\@results); die;
		
		# Clean up 
		system "rm $result_file";
        
        my $i = '0';
        my @closest_blast_matches;
         do {
            
            # code block
            my $hit_ref = $results[$i];
            #$devtools->print_hash($hit_ref); die;
            my $accession = $hit_ref->{'scaffold'};
            print "\n\t Match $accession";
             push (@closest_blast_matches, $accession);
        	$accessions{$accession} = 1;
          
            $i++;
            
            
        } until($i eq 10);
        
        $full_record{$seq_id} = \@closest_blast_matches;
        
        #die;
        # Write the top 10 results file for this input sequence
        #my $file_name = 'output/tabular/' . $seq_id . '-top500_matches.tsv';
        #$fileio->write_file($file_name, \@closest_blast_matches); 

	}

    # Write file of accessions not found in tabular data
	my @missing_accessions = keys %not_found;
	my $missing_accessions_file = 'output/tabular/missing_accessions.tsv';
    $self->{fileio}->write_file($missing_accessions_file, \@missing_accessions); 
	die;

}

1;  # End of FluIO module

