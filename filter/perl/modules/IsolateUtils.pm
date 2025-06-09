package IsolateUtils;

###############################################################################
# Module:       IsolateUtils.pm
# Description:  Provides utility functions for processing influenza isolate
#               records including grouping, consistency checks, and compression.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(order_by_isolate check_consistency compress_isolates);

# ------------------ Subroutine Definitions ------------------

sub order_by_isolate {

    my ($datafile, $ordered_by_isolate_ref, $null_strain_entries_ref) = @_;

    # Read in the file with GenBank entry data
    my @data_file;
    unless ($fileio->read_file($datafile, \@data_file)) {
        die "\n\t # Error reading file '$datafile'\n\n";
    }
    
    # Set up column headers
    my $header_row = shift @data_file;
    chomp $header_row;
    my @header_columns = split("\t", $header_row);
    my %column_indices = map { $header_columns[$_] => $_ } 0 .. $#header_columns;

    # Process each line   
    foreach my $line (@data_file) {
    
        chomp $line;
        my %result;
        @result{@header_columns} = split("\t", $line);

        my $sequenceID  = $result{'sequenceID'};
        my $isolate     = $result{'isolate'};
        my $rec_segment = $result{'rec_segment'};

        # Handle entries with no strain designation
        if (!$isolate || $isolate eq '-') {
            $null_strain_entries_ref->{$sequenceID} = \%result;
            next;
        }

        # Initialize strain entry if it doesn't exist
        my $isolate_ref = $ordered_by_isolate_ref->{$isolate} ||= { 'seq_count' => 0 };

        $isolate_ref->{'seq_count'}++;

        # Initialize segment entry if it doesn't exist
        my $segment_ref = $isolate_ref->{$rec_segment} ||= {};

        # Record this sequence entry
        $segment_ref->{$sequenceID} = \%result;
    }
    
}

sub check_consistency {

    my ($isolate_data, $consistent_ref, $inconsistent_ref) = @_;

    foreach my $isolate_id (keys %$isolate_data) {
    
        my $segments = $isolate_data->{$isolate_id};

        # Variables to store the consistent values for each field
        my %consistent_values = (
            gb_host    => undef,
            gb_country => undef,
            gb_year    => undef,
            gb_month   => undef,
            gb_day     => undef,
        );

        my $is_consistent = 1;

        foreach my $segment_num (@segments) {
        
            my $entries = $segments->{$segment_num};
            
            foreach my $sequenceID (keys %$entries) {
            
                my $entry = $entries->{$sequenceID};

                foreach my $field (keys %consistent_values) {
                    if (defined $consistent_values{$field}) {
                        if ($entry->{$field} ne '-' && $entry->{$field} ne $consistent_values{$field}) {
                            $is_consistent = 0;
                            last;
                        }
                    }
                    else {
                        if ($entry->{$field} ne '-') {
                            $consistent_values{$field} = $entry->{$field};
                        }
                    }
                }

                last unless $is_consistent;
            }

            last unless $is_consistent;
        }

        if ($is_consistent) {
            # Fill in missing values
            foreach my $segment_num (@segments) {
            
                my $entries = $segments->{$segment_num};
                foreach my $sequenceID (keys %$entries) {
                
                    my $entry = $entries->{$sequenceID};
                    foreach my $field (keys %consistent_values) {
                        if ($entry->{$field} eq '-') {
                            $entry->{$field} = $consistent_values{$field};
                        }
                    }
                }
            }
            $consistent_ref->{$isolate_id} = $segments;
        }
        else {
            $inconsistent_ref->{$isolate_id} = $segments;
        }
    }
}

sub compress_isolates {

    my ($isolates_ref, $compressed_ref) = @_;

    foreach my $isolate_id (keys %$isolates_ref) {
		
		#print "\n\t Print ISOLATE $isolate_id";
        
        my $uncompressed_isolate_ref = $isolates_ref->{$isolate_id};
        #$devtools->print_hash($uncompressed_isolate_ref);
        
        my %compressed_isolate;
		my $initialised = undef;
        foreach my $segment_num (@segments) {
        
            if ($uncompressed_isolate_ref->{$segment_num}) {
            
				my $entries = $uncompressed_isolate_ref->{$segment_num};
				#$devtools->print_hash($entries);
				my $longest_entry;
				my $max_length = 0;
				foreach my $sequenceID (keys %$entries) {
					my $entry_ref = $entries->{$sequenceID};
					#$devtools->print_hash($entry_ref);
					
					unless ($initialised) {           		
            			#%compressed_isolate = ;
            			$compressed_isolate{isolate}     = $entry_ref->{isolate};
            			$compressed_isolate{rec_subtype} = $entry_ref->{rec_subtype};
            			$compressed_isolate{gb_subtype}  = $entry_ref->{gb_subtype};
            			$compressed_isolate{rec_subtype} = $entry_ref->{rec_subtype};
            			$compressed_isolate{iso_host}    = $entry_ref->{iso_host};
            			$compressed_isolate{pubmed_id}   = $entry_ref->{pubmed_id};
            			$compressed_isolate{iso_day}     = $entry_ref->{iso_day};
            			$compressed_isolate{iso_month}   = $entry_ref->{iso_month};
            			$compressed_isolate{iso_year}    = $entry_ref->{iso_year};
            			$compressed_isolate{iso_country} = $entry_ref->{iso_country};
             			$compressed_isolate{iso_source}  = $entry_ref->{iso_source};
            			$compressed_isolate{lab_host}    = $entry_ref->{lab_host};
            			$compressed_isolate{'length'}    = $entry_ref->{'length'};
            			$compressed_isolate{rec_segment} = $entry_ref->{rec_segment};
            			$initialised = 1;
            		}
					
					if ($entry_ref->{length} > $max_length) {
						$max_length = $entry_ref->{length};
						$longest_entry = $entry_ref;
					}
				}

				# Store the longest entry for this segment
				my $key = 'segment' . $segment_num . '_accession';
				#$devtools->print_hash($longest_entry); die;
				$compressed_isolate{$key} = $longest_entry->{sequenceID};

            }
        }
        
        $compressed_ref->{$isolate_id} = \%compressed_isolate;
    
    }
}

1;  # End of IsolateUtils module
