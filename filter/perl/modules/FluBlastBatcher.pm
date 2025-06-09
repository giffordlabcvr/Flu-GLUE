package FluBlastBatcher;

###############################################################################
# Module:       FluBlastBatcher.pm
# Description:  Provides functionality for BLASTing influenza genome segments
#               against a local RSL library and recording top hits.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(batch_blast_against_rsl);

# ------------------ Subroutine Definitions ------------------

sub batch_blast_against_rsl {
    my ($infile, $seqio, $fileio, $blast_obj, $lib_path, $tmp_path, $output_path) = @_;

	my ($infile) = @_;

	unless ($lib_path and $tmp_path and $output_path) { die; }

	# Read FASTA input
	my @input_fasta;
	$seqio->read_fasta($infile, \@input_fasta);

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
		$fileio->write_text_to_file($query_file, $fasta);
		my $result_file = $tmp_path . "/$seq_id.blast_result.tsv";

		# Execute the call to BLAST and parse the results
		$blast_obj->blast('/Users/rob/blast/bin/blastn', $lib_path, $query_file, $result_file);
		my @results;
		$blast_obj->parse_tab_format_results($result_file, \@results);
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
    $fileio->write_file($missing_accessions_file, \@missing_accessions); 
	die;

}

1;  # End of FluBlastBatcher module
