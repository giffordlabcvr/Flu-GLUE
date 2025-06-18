#!/usr/bin/perl -w
############################################################################
# Module:      FluIsolateUtils.pm 
# Description: 
# History:     Rob Gifford January 2024: Creation
############################################################################
package FluIsolateUtils;

############################################################################
# Import statements
############################################################################
use strict;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new FluIsolateUtils.pm 'object'
#***************************************************************************
sub new {

    my ($class, $fileio, $num_segments, $species, $devtools) = @_;

    my $self = {
        fileio        => $fileio,
        devtools      => $devtools,
        num_segments  => $num_segments,
        species       => $species,
    };
    bless $self, $class;
    return $self;
}


############################################################################
# Subroutine Definitions
############################################################################

#***************************************************************************
# Subroutine:  order_by_isolate
# Description: Index influenza GenBank entries by strain name
# Arguments: $datafile (array containing GenBank entry table)
#            $ordered_by_isolate_ref (hash to store GenBank entries indexed by strain)
#***************************************************************************
sub order_by_isolate {

    my ($self, $datafile, $ordered_by_isolate_ref, $null_strain_entries_ref) = @_;

    # Read in the file with GenBank entry data
    my @data_file;
    unless ($self->{fileio}->read_file($datafile, \@data_file)) {
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

#***************************************************************************
# Subroutine:  check_consistency
# Description: check the consistency of isolate-associated data across influenza segments
#***************************************************************************
sub check_consistency {

    my ($self, $isolate_data, $consistent_ref, $inconsistent_ref) = @_;

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
 
		for (my $segment_num = 1; $segment_num <= $self->{num_segments}; $segment_num++) {

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
			for (my $segment_num = 1; $segment_num <= $self->{num_segments}; $segment_num++) {

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

#***************************************************************************
# Subroutine:  compress_isolates
# Description: For each isolate, select the longest sequence entry for each segment
# Arguments: $isolates_ref (hash reference containing consistent isolates)
#            $compressed_ref (hash reference to store compressed isolates)
#***************************************************************************
sub compress_isolates {
    my ($self, $isolates_ref, $compressed_ref) = @_;

    my $species = $self->{species};
    my $devtools = $self->{devtools};

    # Define which segments are used for serotype determination
    my %serotype_segment = (
        iav => [4, 6],
        ibv => [4],
        icv => [4],
        idv => [4],
    );

    foreach my $isolate_id (keys %$isolates_ref) {
    
        my $uncompressed_isolate_ref = $isolates_ref->{$isolate_id};
        my %compressed_isolate;
        my $initialised = undef;
        my @lineages;
	    my %serotype_data;

        for (my $segment_num = 1; $segment_num <= $self->{num_segments}; $segment_num++) {
           
            if ($uncompressed_isolate_ref->{$segment_num}) {
            
                my $entries = $uncompressed_isolate_ref->{$segment_num};
                my $longest_entry;
                my $max_length = 0;


				foreach my $sequenceID (keys %$entries) {
					my $entry_ref = $entries->{$sequenceID};
					
					#$devtools->print_hash($entry_ref); die;

					unless ($initialised) {
						$compressed_isolate{isolate}        = $entry_ref->{isolate};
						$compressed_isolate{iso_host}       = $entry_ref->{gb_host};
						$compressed_isolate{pubmed_id}      = $entry_ref->{pubmed_id};
						$compressed_isolate{iso_day}        = $entry_ref->{gb_day};
						$compressed_isolate{iso_month}      = $entry_ref->{gb_month};
						$compressed_isolate{iso_year}       = $entry_ref->{gb_year};
						$compressed_isolate{iso_country}    = $entry_ref->{gb_country};
						$compressed_isolate{iso_place_name} = $entry_ref->{gb_place_name};
						$compressed_isolate{lab_host}       = $entry_ref->{gb_lab_host};
						$compressed_isolate{sample_type}    = $entry_ref->{gb_source};

						# Derive origin_type and host
						my $host_str = $entry_ref->{gb_host} // '';
						my $isolate_str = $entry_ref->{isolate} // '';

						if (
							$host_str =~ /^environment$/i ||
							$host_str =~ /environmental|air|water|surface/i ||
							$isolate_str =~ m#^[A-Z]/environment/#i
						) {
							$compressed_isolate{origin_type} = 'environment';
							$compressed_isolate{host}        = '-';
						} else {
							$compressed_isolate{origin_type} = 'host';
							$compressed_isolate{host}        = $host_str ne '' ? $host_str : '-';
						}

						$initialised = 1;
					}


					my $prefer_this = 0;

					# Prefer entries with serotype data if this segment is designated for serotyping
					if (grep { $_ == $segment_num } @{ $serotype_segment{$species} || [] }) {
						if (defined $entry_ref->{rec_subtype} && $entry_ref->{rec_subtype} ne '-') {
							$prefer_this = 1;
						}
					}

					if (
						!$longest_entry ||
						$prefer_this ||
						($entry_ref->{length} > $max_length && !$prefer_this)
					) {
						$max_length = $entry_ref->{length};
						$longest_entry = $entry_ref;
					}
					
					
				}



                # Assign accession number for the segment
                my $key = 'segment' . $segment_num . '_accession';
                $compressed_isolate{$key} = $longest_entry->{sequenceID};

                # Capture subtype info if segment is used for serotyping
                if (grep { $_ == $segment_num } @{ $serotype_segment{$species} || [] }) {
                    $serotype_data{rec_serotype}  ||= $longest_entry->{rec_subtype}  if $longest_entry->{rec_subtype} ne '-';
                    $serotype_data{gb_serotype}   ||= $longest_entry->{gb_subtype}   if $longest_entry->{gb_subtype}  ne '-';
                    $serotype_data{mlca_serotype} ||= $longest_entry->{rec_subtype} if $longest_entry->{rec_subtype} ne '-';
                }

                # Accumulate for genome_lineage (may contain dashes)
                push @lineages, ($longest_entry->{rec_subtype} || '-');
            } else {
                push @lineages, '-';  # No segment, placeholder
            }
        }
        
        
        

        $compressed_isolate{gb_serotype}    = $serotype_data{gb_subtype}   || '-';
        $compressed_isolate{rec_serotype}   = $serotype_data{rec_subtype}  || '-';
        $compressed_isolate{mlca_serotype}  = $serotype_data{mlca_serotype}|| '-';

        # Create genome_lineage from all segments' recogniser-based subtypes
        my $lineage_string = join('|', @lineages);
        $compressed_isolate{genome_lineage} = $lineage_string;

		#$devtools->print_hash(\%compressed_isolate);

        $compressed_ref->{$isolate_id} = \%compressed_isolate;
    }
}




#***************************************************************************
# Subroutine:  check_completeness
# Description: sort isolates into those that are complete versus incomplete
#***************************************************************************
sub check_completeness {

	my ($self, $ordered_by_isolate_ref, $complete_ref, $incomplete_ref, $species) = @_;
		
	# Iterate through isolates
	my @isolate_ids = keys %$ordered_by_isolate_ref;
	foreach my $isolate_id (@isolate_ids) {
	
		my $isolate_ref = $ordered_by_isolate_ref->{$isolate_id};
		#$devtools->print_hash($isolate_ref); die;
	
	    my %segment_presence;
		for (my $segment_num = 1; $segment_num <= $self->{num_segments}; $segment_num++) {
	    	if ($isolate_ref->{$segment_num}) {	    		
	    		$segment_presence{$segment_num} = 1;
	    	}
	    }
	    
	    my $seg_coverage = scalar keys %segment_presence;
	    
	    # If its IAV or IBV check if eight segments present
	    if ($species eq 'iav' or $species eq 'ibv') {
	    	
	    	if ($seg_coverage eq 8) {
	    		$complete_ref->{$isolate_id} = $isolate_ref;
	    	}
	    	else {
	    		$incomplete_ref->{$isolate_id} = $isolate_ref;	    	
	    	}
	    }
	    # If its ICV or IDV check if seven segments present
	    if ($species eq 'icv' or $species eq 'idv') {
	    
	    	if ($seg_coverage eq 7) {
	    		$complete_ref->{$isolate_id} = $isolate_ref;
	    	}
	    	else {
	    		$incomplete_ref->{$isolate_id} = $isolate_ref;	    	
	    	}
	    	
	    }
	    
	}
	
}

1;  # End of FluIsolateUtils module

