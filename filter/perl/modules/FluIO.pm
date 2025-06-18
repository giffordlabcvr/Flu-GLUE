#!/usr/bin/perl -w
############################################################################
# Module:      FluIO.pm 
# Description: 
# History:     Rob Gifford January 2024: Creation
############################################################################
package FluIO;

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

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new FluIO.pm 'object'
#***************************************************************************
sub new {

    my ($class, $fileio, $num_segments, $module_path, $glue_out_dir, $devtools) = @_;

	# Member variables
    my $self = {
        fileio    => $fileio,
        devtools  => $devtools,
        num_segments  => $num_segments,
        module_path   => $module_path,
        glue_out_dir  => $glue_out_dir,

    };
    bless $self, $class;
    return $self;
}

############################################################################
# Subroutine Definitions
############################################################################

#***************************************************************************
# Subroutine:  write_isolate_data_table
# Description: Write isolate data to a tabular file
# Arguments: $isolates_ref (hash reference containing uncompressed isolates data:
#            (i.e. with redundant segment sequences included)
#            $output_file (output file path)
#            $species (species type to determine number of segments)
#***************************************************************************
sub write_isolate_data_table {

    my ($self, $isolates_ref, $output_file, $species) = @_;

    # Create header line
    my @header_fields = qw [ isolate gb_serotype rec_serotype mlca_serotype genome_lineage
                             iso_year iso_month iso_day
                             origin_type iso_host sample_type
                             iso_country iso_place_name
                             segment1_accession segment2_accession segment3_accession segment4_accession 
                             segment5_accession segment6_accession segment7_accession ];

	my $devtools = $self->{devtools}; 


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
            gb_serotype => '-',
            rec_serotype => '-',
            mlca_serotype => '-',
            genome_lineage => '-',
            iso_year => '-',
            iso_month => '-',
            iso_day => '-',
            origin_type => '-',
            iso_host => '-',
            sample_type => '-',
            iso_country => '-',
            iso_place_name => '-',
            
            segment1_accession => '-',
            segment2_accession => '-',
            segment3_accession => '-',
            segment4_accession => '-',
            segment5_accession => '-',
            segment6_accession => '-',
            segment7_accession => '-',
            segment8_accession => '-'
        );
		
		#$devtools->print_hash($isolate_ref);	die;
		
		# Fill in the data fields
		# Note: these values will be written over with each iteration but should in any case be identical across segments
		$row_data{iso_year} = $isolate_ref->{iso_year} if $isolate_ref->{iso_year};
		$row_data{iso_month} = $isolate_ref->{iso_month} if $isolate_ref->{iso_month};
		$row_data{iso_day} = $isolate_ref->{iso_day} if $isolate_ref->{iso_day};
		$row_data{origin_type} = $isolate_ref->{origin_type} if $isolate_ref->{origin_type};
		$row_data{iso_host} = $isolate_ref->{iso_host} if $isolate_ref->{iso_host};
		$row_data{sample_type} = $isolate_ref->{sample_type} if $isolate_ref->{sample_type};
		$row_data{iso_country} = $isolate_ref->{iso_country} if $isolate_ref->{iso_country};
		$row_data{iso_place_name} = $isolate_ref->{iso_place_name} if $isolate_ref->{iso_place_name};
		$row_data{gb_serotype} = $isolate_ref->{gb_serotype} if $isolate_ref->{gb_serotype};
		$row_data{rec_serotype} = $isolate_ref->{rec_serotype} if $isolate_ref->{rec_serotype};
		$row_data{mlca_serotype} = $isolate_ref->{mlca_serotype} if $isolate_ref->{mlca_serotype};
		$row_data{genome_lineage} = $isolate_ref->{genome_lineage} if $isolate_ref->{genome_lineage};
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

    $self->{fileio}->write_file($output_file, \@output);
}


#***************************************************************************
# Subroutine:  export_gb_entries_from_isolate_hash
# Description: export a tabular file of GenBank entry details from isolate hash
#***************************************************************************
sub export_gb_entries_from_isolate_hash {

    my ($self, $output_file, $isolates_ref) = @_;

    # Create header line
    my @output;
	my $devtools = $self->{devtools};
	
    my @header_fields = qw [ sequenceID isolate gb_segment rec_segment 
                             gb_serotype rec_serotype                               
                             iso_year iso_month iso_day
                             iso_host iso_country iso_place_name
                             pubmed_id length ];
                             
    my $header = join("\t", @header_fields);
    push (@output, "$header\n");

    foreach my $isolate_id (keys %$isolates_ref) {
    
		my $isolate_details_ref = $isolates_ref->{$isolate_id};
		for (my $segment_num = 1; $segment_num <= $self->{num_segments}; $segment_num++) {
		
		    my $segment_seqs_ref = $isolate_details_ref->{$segment_num};
		    if ($segment_seqs_ref) {
		    
				my @sequence_ids = keys %$segment_seqs_ref;
				
				foreach my $sequence_id (@sequence_ids) {
			
					# Print the row data to the file
					my $gb_entry_details = $segment_seqs_ref->{$sequence_id};
					$devtools->print_hash($gb_entry_details);
					
					my @row_values;
					push (@row_values, $sequence_id);
					my $isolate = '-';
					my $gb_segment = '-';
					my $rec_segment = '-';
					my $gb_serotype = '-';
					my $rec_serotype = '-';
					my $iso_year = '-';
					my $iso_month = '-';
					my $iso_day = '-';
					my $origin_type = '-';
					my $iso_host = '-';
					my $sample_type = '-';
					my $iso_country = '-';
					my $iso_place_name = '-';
					my $pubmed_id = '-';
					my $length = $gb_entry_details->{length};

					unless ($length) { die "\n\t Sequence length is not defined\n\n"; } # Every sequence should have a length

					if ($gb_entry_details->{isolate}) {
						$gb_segment =  $gb_entry_details->{isolate};		
					}
					if ($gb_entry_details->{gb_segment}) {
						$gb_segment =  $gb_entry_details->{gb_segment};		
					}
					if ($gb_entry_details->{rec_segment}) {
						$rec_segment =  $gb_entry_details->{rec_segment};		
					}										
					if ($gb_entry_details->{gb_subtype}) {
						$gb_serotype =  $gb_entry_details->{gb_subtype};		
					}
					if ($gb_entry_details->{rec_subtype}) {
						$rec_serotype =  $gb_entry_details->{rec_subtype};		
					}
					if ($gb_entry_details->{gb_year}) {
						$iso_year =  $gb_entry_details->{gb_year};		
					}
					if ($gb_entry_details->{gb_month}) {
						$iso_month =  $gb_entry_details->{gb_month};		
					}
					if ($gb_entry_details->{gb_day}) {
						$iso_day =  $gb_entry_details->{gb_day};		
					}
					if ($gb_entry_details->{gb_host}) {
						$iso_host =  $gb_entry_details->{'gb_host'};		
					}
					if ($gb_entry_details->{gb_country}) {
						$iso_country =  $gb_entry_details->{'gb_country'};		
					}
					if ($gb_entry_details->{gb_place_name}) {
						$iso_place_name =  $gb_entry_details->{'gb_place_name'};		
					}
					if ($gb_entry_details->{pubmed_id}) {
						$pubmed_id =  $gb_entry_details->{'pubmed_id'};		
					}
					if ($gb_entry_details->{'length'}) {
						$length =  $gb_entry_details->{'length'};		
					}

					push (@row_values, $isolate_id);
					push (@row_values, $gb_segment);
					push (@row_values, $rec_segment);
					push (@row_values, $gb_serotype);
					push (@row_values, $rec_serotype);
					push (@row_values, $iso_year);
					push (@row_values, $iso_month);
					push (@row_values, $iso_day);
					push (@row_values, $iso_host);
					push (@row_values, $iso_country);
					push (@row_values, $iso_place_name);
					push (@row_values, $pubmed_id);
					push (@row_values, $length);
					
					my $row = join("\t", @row_values);
					push (@output, "$row\n");
				
				}
			
			}

		}
				
    }
    
	$self->{fileio}->write_file($output_file, \@output);
}

#***************************************************************************
# Subroutine:   create_glue_download_program_logic
# Description:  Generates GLUE-compatible module and command logic for
#               downloading curated isolate sequences from NCBI.
#***************************************************************************
sub create_glue_download_program_logic {

	my ($self, $isolates_ref, $species, $descriptorG, $descriptorM) = @_;    
	
	
	# GLUE FILE STRINGS
	my @delete_text;
	my @create_text;	
	my @import_text;	
	my @export_text;	
	my $glue_file_name = $species . 'Download' . $descriptorM . '.glue';
	my $glue_file_path = $self->{glue_out_dir} . '/' . $glue_file_name;

	my @isolate_ids = keys %$isolates_ref;
	my $isolate_count = scalar @isolate_ids;

	# Iterate through segments at the top because we download each separately
	for (my $segment_num = 1; $segment_num <= $self->{num_segments}; $segment_num++) {

		# Set up the source name and import line for GLUE file
		my @import_source;	
		my @download_source;	
		my $source_name = $species . "-" . $descriptorG . "-segment" . $segment_num;
		my $import_line = "\n\t#import source sources/$source_name";
		my $download_line = "\n\t#run file $glue_file_path ";
		push (@import_source, $import_line);	
		push (@download_source, $download_line);

		# Set up a module text string
		my $module_name = $species . 'NcbiImporter' . $descriptorM . "Segment" . $segment_num;
		my @module_text;
		my $header = "<!-- NCBI import module for $species '$descriptorG' segment $segment_num source -->\n<ncbiImporter>\n";
		$header .= "\t<sequenceFormat>GENBANK_XML</sequenceFormat>\n\t<sourceName>$source_name</sourceName>\n\t<specificPrimaryAccessions>";		
		my $tail = "\n\t</specificPrimaryAccessions>\n\t<sequenceIdField>PRIMARY_ACCESSION</sequenceIdField>\n</ncbiImporter>";
		print "\n\t CREATING DOWNLOAD MODULES FOR '$descriptorG' SEGMENT '$segment_num':";
		push (@module_text, $header);

		# Iterate through isolates
		my $added = 0;
		foreach my $isolate_id (@isolate_ids) {
	
			# Get isolate details
			my $isolate_ref = $isolates_ref->{$isolate_id};
			#$devtools->print_hash($isolate_ref); die;
			my $key = 'segment' . $segment_num . '_accession';
			if ($isolate_ref->{$key}) {			
				my $accession = $isolate_ref->{$key};
				my $line = "\n\t\t<primaryAccession>$accession</primaryAccession> <!-- $isolate_id -->";
				push (@module_text, $line);
				$added++;
			}

		}
		
		if ($added) {
			push (@module_text, $tail);
			my $file_name = $module_name . '.xml';
			$self->{fileio}->write_file($file_name, \@module_text);

			my $delete_line = "\n\tdelete module $module_name"; 
			push (@delete_text, $delete_line);
			my $create_line = "\n\tcreate module -f $self->{module_path}/$species/$module_name" . '.xml'; 
			push (@create_text, $create_line);

			my $import_cmd_line = "\n\tmodule $module_name import";
			push (@import_text, $import_cmd_line);

			my $export_cmd_line = "\n\texport source $source_name";
			push (@export_text, $export_cmd_line);
		
		}
	}

    # ASSEMBLE GLUE FILE and write as output
    my $uc_species = uc $species;
    my $header_comment = "\t# This GLUE script downloads from GenBank $uc_species sequences representing $isolate_count isolates\n\n";

    my $delete_comment = "\t# Delete previous modules if they exist";
    my $create_comment = "\t# Create the download modules";

	my @set_up_text = ( "$delete_comment@delete_text\n\n", "$create_comment@create_text\n\n" );	
	unshift (@set_up_text, $header_comment);

    my $import_comment .= "\t# Download segment sequences - each is downloaded to a separate source";
    my $export_comment .= "\t# Export segment sequences";
	my @final_output = ( "@set_up_text", "$import_comment@import_text\n\n", "$export_comment@export_text\n\n" );
	
	# Write to current directory
	$self->{fileio}->write_file($glue_file_name, \@final_output);

}


1;  # End of FluIO module

