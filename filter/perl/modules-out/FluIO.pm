package FluIO;

###############################################################################
# Module:       FluIO.pm
# Description:  Provides input/output utilities for exporting influenza isolate
#               data to tabular and GenBank-compatible formats.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(write_isolate_data_table export_gb_entries_from_isolate_hash);

# ------------------ Subroutine Definitions ------------------

sub write_isolate_data_table {

    my ($isolates_ref, $output_file, $species) = @_;

    # Create header line
    my @header_fields = qw [ isolate gb_subtype rec_subtype  
                             iso_year iso_month iso_day
                             iso_host iso_source iso_country 
                             segment1_accession segment2_accession segment3_accession segment4_accession 
                             segment5_accession segment6_accession segment7_accession ];

    if ($species eq 'iav' or $species eq 'ibv') {
        push(@header_fields, 'segment8_accession');
    }
    my @output;
    my $header = join("\t", @header_fields);
    push (@output, "$header\n");

    # Iterate through isolates
    foreach my $isolate_id (keys %$isolates_ref) {
    
        my $isolate_ref = $isolates_ref->{$isolate_id};

        # Initialize a hash to store row data
        my %row_data = (
            isolate => $isolate_id,
            gb_subtype => '-',
            rec_subtype => '-',
            iso_year => '-',
            iso_month => '-',
            iso_day => '-',
            iso_host => '-',
            iso_source => '-',
            iso_country => '-',
            segment1_accession => '-',
            segment2_accession => '-',
            segment3_accession => '-',
            segment4_accession => '-',
            segment5_accession => '-',
            segment6_accession => '-',
            segment7_accession => '-',
            segment8_accession => '-'
        );
		
		#$devtools->print_hash($isolate_ref);	
		
		# Fill in the data fields
		# Note: these values will be written over with each iteration but should in any case be identical across segments
		$row_data{iso_year} = $isolate_ref->{iso_year} if $isolate_ref->{iso_year};
		$row_data{iso_month} = $isolate_ref->{iso_month} if $isolate_ref->{iso_month};
		$row_data{iso_day} = $isolate_ref->{iso_day} if $isolate_ref->{iso_day};
		$row_data{iso_host} = $isolate_ref->{iso_host} if $isolate_ref->{iso_host};
		$row_data{iso_source} = $isolate_ref->{iso_source} if $isolate_ref->{iso_source};
		$row_data{iso_country} = $isolate_ref->{iso_country} if $isolate_ref->{iso_country};
		$row_data{gb_subtype} = $isolate_ref->{gb_subtype} if $isolate_ref->{gb_subtype};
		$row_data{rec_subtype} = $isolate_ref->{rec_subtype} if $isolate_ref->{rec_subtype};
		$row_data{segment1_accession} = $isolate_ref->{segment1_accession} if $isolate_ref->{segment1_accession};
		$row_data{segment2_accession} = $isolate_ref->{segment2_accession} if $isolate_ref->{segment2_accession};
		$row_data{segment3_accession} = $isolate_ref->{segment3_accession} if $isolate_ref->{segment3_accession};
		$row_data{segment4_accession} = $isolate_ref->{segment4_accession} if $isolate_ref->{segment4_accession};
		$row_data{segment5_accession} = $isolate_ref->{segment5_accession} if $isolate_ref->{segment5_accession};
		$row_data{segment6_accession} = $isolate_ref->{segment6_accession} if $isolate_ref->{segment6_accession};
		$row_data{segment7_accession} = $isolate_ref->{segment7_accession} if $isolate_ref->{segment7_accession};
		$row_data{segment8_accession} = $isolate_ref->{segment8_accession} if $isolate_ref->{segment8_accession};

        # Print the row data to the file
        my @row_values;
        foreach my $field (@header_fields) {
        	
        	my $value = $row_data{$field};
        	unless ($value) {
        		print "\n\t Got null value for field $field\n\n\n"; die;
        	}
        	push (@row_values, $value);
        	
        }
        my $row = join("\t", @row_values);
        push (@output, "$row\n");
    }

    $fileio->write_file($output_file, \@output);
}

sub export_gb_entries_from_isolate_hash {

    my ($output_file, $isolates_ref) = @_;

    # Create header line
    my @output;
    my @header_fields = qw [ sequenceID isolate gb_segment rec_segment 
                             gb_subtype rec_subtype                               
                             iso_year iso_month iso_day
                             iso_host iso_source iso_country iso_place_name
                             pubmed_id length ];
    my $header = join("\t", @header_fields);
    push (@output, "$header\n");

    foreach my $isolate_id (keys %$isolates_ref) {
    
		my $isolate_details_ref = $isolates_ref->{$isolate_id};
		foreach my $segment_num (@segments) {
		
		    my $segment_seqs_ref = $isolate_details_ref->{$segment_num};
		    if ($segment_seqs_ref) {
		    
				my @sequence_ids = keys %$segment_seqs_ref;
				
				foreach my $sequence_id (@sequence_ids) {
			
					# Print the row data to the file
					my $gb_entry_details = $segment_seqs_ref->{$sequence_id};
					#$devtools->print_hash($gb_entry_details); die;
					
					my @row_values;
					push (@row_values, $sequence_id);
					my $gb_segment = '-';
					my $rec_segment = '-';
					my $gb_subtype = '-';
					my $rec_subtype = '-';
					my $iso_year = '-';
					my $iso_month = '-';
					my $iso_day = '-';
					my $iso_host = '-';
					my $iso_source = '-';
					my $iso_country = '-';
					my $iso_place_name = '-';
					my $pubmed_id = '-';

					if ($gb_entry_details->{gb_segment}) {
						$gb_segment =  $gb_entry_details->{gb_segment};		
					}
					if ($gb_entry_details->{rec_segment}) {
						$rec_segment =  $gb_entry_details->{rec_segment};		
					}										
					if ($gb_entry_details->{gb_subtype}) {
						$gb_subtype =  $gb_entry_details->{gb_subtype};		
					}
					if ($gb_entry_details->{rec_subtype}) {
						$rec_subtype =  $gb_entry_details->{rec_subtype};		
					}
					if ($gb_entry_details->{iso_year}) {
						$iso_year =  $gb_entry_details->{iso_year};		
					}
					if ($gb_entry_details->{iso_month}) {
						$iso_month =  $gb_entry_details->{iso_month};		
					}
					if ($gb_entry_details->{iso_day}) {
						$iso_day =  $gb_entry_details->{iso_day};		
					}
					if ($gb_entry_details->{iso_host}) {
						$iso_host =  $gb_entry_details->{iso_host};		
					}
					if ($gb_entry_details->{iso_source}) {
						$iso_source =  $gb_entry_details->{iso_source};		
					}
					if ($gb_entry_details->{iso_country}) {
						$iso_country =  $gb_entry_details->{iso_country};		
					}
					if ($gb_entry_details->{iso_place_name}) {
						$iso_place_name =  $gb_entry_details->{iso_place_name};		
					}
					if ($gb_entry_details->{pubmed_id}) {
						$pubmed_id =  $gb_entry_details->{pubmed_id};		
					}
					
					push (@row_values, $sequence_id);
					push (@row_values, $isolate_id);
					push (@row_values, $gb_segment);
					push (@row_values, $rec_segment);
					push (@row_values, $gb_subtype);
					push (@row_values, $rec_subtype);
					push (@row_values, $iso_year);
					push (@row_values, $iso_month);
					push (@row_values, $iso_day);
					push (@row_values, $iso_source);
					push (@row_values, $iso_country);
					push (@row_values, $iso_place_name);
					push (@row_values, $pubmed_id);
					push (@row_values, $gb_entry_details->{length});
					
					my $row = join("\t", @row_values);
					push (@output, "$row\n");
				
				}
			
			}

		}
				
    }
    
	$fileio->write_file($output_file, \@output);
}

1;  # End of FluIO module
