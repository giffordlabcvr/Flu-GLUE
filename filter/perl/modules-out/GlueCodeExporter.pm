package GlueCodeExporter;

###############################################################################
# Module:       GlueCodeExporter.pm
# Description:  Generates GLUE-compatible module and command logic for
#               downloading curated isolate sequences from NCBI.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(create_glue_download_program_logic);

# ------------------ Subroutine Definitions ------------------

sub create_glue_download_program_logic {

	my ($isolates_ref, $species, $descriptorG, $descriptorM) = @_;    
	
	
	# GLUE FILE STRINGS
	my @delete_text;
	my @create_text;	
	my @import_text;	
	my @export_text;	
	my $glue_file_name = $species . 'Download' . $descriptorM . '.glue';
	my $glue_file_path = $glue_dl_code_directory . '/' . $glue_file_name;

	my @isolate_ids = keys %$isolates_ref;
	my $isolate_count = scalar @isolate_ids;

	# Iterate through segments at the top because we download each separately
	foreach my $segment_num (@segments) {

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
			$fileio->write_file($file_name, \@module_text);

			my $delete_line = "\n\tdelete module $module_name"; 
			push (@delete_text, $delete_line);
			my $create_line = "\n\tcreate module -f $module_path/$species/$module_name" . '.xml'; 
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
	$fileio->write_file($glue_file_name, \@final_output);

}

1;  # End of GlueCodeExporter module
