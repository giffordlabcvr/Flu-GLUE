#!/usr/bin/perl -w
############################################################################
# Script:      flu-gb-filter-tool.pl 
# Description: a collection of tools for working with sequences + data
# History:     Version 1.0 Creation: Rob J Gifford 2024
############################################################################

# set up library of local modules for this program (base functions)
use lib './modules/'; 

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use Getopt::Long;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base modules
use Console;
use DevTools;
use FileIO;
use SeqIO; 

# Interface modules
use BLAST;

############################################################################
# Globals
############################################################################

# Version number
my $program_version = '0.1 beta';

# Segment array (for ICV and IDV the '8' is removed via 'pop')
my @segments = qw [ 1 2 3 4 5 6 7 8 ];

# Paths for output (these should reflect current organisation of Flu-GLUE)
my $glue_dl_code_directory = 'glue/download/curated';
my $module_path = "modules/download/curated/";

# Input file paths
my $country_region_data_path = '../tabular/iso_regions/all.tsv';
my $lib_path = '/Users/rob/Sources/flu/blastdb/iav/rsl-flu-epi.fna'; 
my $tmp_path = './tmp/';
my $output_path = './';

############################################################################
# Instantiations for program 'classes' (PERL's Object-Oriented Emulation)
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();
my $seqio      = SeqIO->new();
my $blast_obj  = BLAST->new();

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 -m=[1] -i=[infile] -s=[species]\n\n\n";

############################################################################
# Main program
############################################################################

# Run script
main();

# Exit program
exit;

############################################################################
# Subroutines
############################################################################

#***************************************************************************
# Subroutine:  main
# Description: top level handler fxn
#***************************************************************************
sub main {

	# Define variables to hold the option values
	my ($help, $version, $mode, $species, $file1, $file2);
	
	# Read in options using GetOpt::Long
	GetOptions(
		'help!'      => \$help,
		'version!'   => \$version,
		'mode|m=i'   => \$mode,
		'species|s=s'=> \$species,
		'infile|i=s' => \$file1,
		'infile2|f=s'=> \$file2,
	) or die "Error in command line arguments";

	if ($help) { # Show help page
		show_help_page();  
	}
	elsif ($version)  { 
		print "\n\t # flu-gb-filter-tool.pl version $program_version\n\n";  
	}
	elsif ($mode)  {
	
		show_title();
	    unless ($file1) {	    
	    	print "\n\t  Please supply an input file for mode '$mode' using the '-i' option\n";
	    	print "\n\t  For more details on how to run, use '-h' (help)\n";
	    	die $USAGE;
	    }
		if (($mode eq 1) and $file1) {		    
		    create_isolate_database($file1, $species);
		}
		elsif (($mode eq 2) and $file1) {
			select_isolate_subsets($file1, $species);
		}
		elsif (($mode eq 3) and $file1) {
            my %ordered_by_isolate;
			batch_blast_against_rsl($file1);
		}
		elsif (($mode eq 4) and $file1) {
			remove_sequences_from_source_directory($file1, $file2);
		}
		elsif (($mode eq 5) and $file1) {
			summarise_influenza_sequence_set($file1);
		}
		elsif (($mode eq 6) and $file1) {
			convert_iav_strain_names_to_data_fields($file1);
		}

	}
	else {
		die $USAGE;
	}
	print "\n\n\t # Exit\n\n\n";
}

############################################################################
# Compile isolate database
############################################################################

#***************************************************************************
# Subroutine:  create_isolate_database
# Description: derive isolate-oriented data from a table of influenzavirus segment sequences
#***************************************************************************
sub create_isolate_database {

	my ($datafile, $species) = @_;

    # Adjust segments based on species
	adjust_segments($species);

    # Order sequences by strain 
	my %ordered_by_isolate;
	my %null_isolate_entries;
	order_by_isolate($datafile, \%ordered_by_isolate, \%null_isolate_entries);
    #$devtools->print_hash(\%ordered_by_isolate); die;

	# Write the null isolate entries to file
	#my $outfile_null_isolates = $species . "_null_isolate-id_entries.tsv";
	#export_gb_entries($outfile_null_isolates, \%null_isolate_entries);

    # Check consistency of isolate-associated, demographic data across different segment sequences
    my %consistent;
    my %inconsistent;
    check_consistency(\%ordered_by_isolate, \%consistent, \%inconsistent);
	my $outfile_inconsistent = $species . "_inconsistent_isolate_entries.tsv";
	# Write the inconsistent entries to file
	export_gb_entries_from_isolate_hash($outfile_inconsistent, \%inconsistent);
	
    # For isolates with consistent data, check completeness
    my %incomplete;
    my %complete;
    check_completeness(\%consistent, \%complete, \%incomplete, $species); 
    #$devtools->print_hash(\%incomplete); die;

    # Compress sequence sets for both incomplete & complete isolates.
    # In other words: where multiple sequences are available for an isolate segment, select one
    my %compressed_complete;
    my %compressed_incomplete;
    compress_isolates(\%complete, \%compressed_complete); 
    compress_isolates(\%incomplete, \%compressed_incomplete);
    #$devtools->print_hash(\%compressed_incomplete); die;

    # TODO: Determine the complete genome subtype for IAV isolates with all eight segments
    if ($species eq 'iav') {
    	my %counts;
    	#determine_iav_genome_lineage(\%compressed_complete, \%counts);
    }

    # Write isolate-ordered tabular data for complete genome sequences
    my $outfile_complete = $species . "_complete_isolates.tsv";
    write_isolate_data_table(\%compressed_complete, $outfile_complete, $species); 
    
    # Write isolate-ordered tabular data for isolates lacking complete genome info
    my $outfile_incomplete = $species . "_incomplete_isolates.tsv";
    write_isolate_data_table(\%compressed_incomplete, $outfile_incomplete, $species); 
	
	# Give option to export GLUE program logic as well as isolate db
	my $num_complete_isolates = scalar keys %compressed_complete;
	if ($num_complete_isolates <= 1000) {
	
		my $code_export_question = "\n\t  Do you want to export GLUE modules and code for downloading these isolates?";
		my $export_response = $console->ask_yes_no_question($code_export_question);
		if ($export_response eq 'y') {
			
			# Complete genome isolates 
			{
				# create GLUE program logic for downloading complete genome isolates
				my $descriptorG = 'ncbi-curated'; # For source name 
				my $descriptorM = 'NcbiCurated';  # For module name
				create_glue_download_program_logic(\%compressed_complete, $species, $descriptorG, $descriptorM);
			
			}
			# Partial genome isolates 
			{
				# create GLUE program logic for downloading complete genome isolates
				my $descriptorG = 'ncbi-curated-partial'; # For source name 
				my $descriptorM = 'NcbiCuratedPartial';  # For module name
				create_glue_download_program_logic(\%compressed_incomplete, $species, $descriptorG, $descriptorM);
			
			}
		
		}
		
	}

}

#***************************************************************************
# Subroutine:  select_isolate_subsets
# Description: select subsets of isolate sequences, stratified by their 
#***************************************************************************
sub select_isolate_subsets {

	my ($datafile, $species) = @_;

    # Adjust segments based on species
	adjust_segments($species);

    # Read in the isolate data
	my @isolates;
	$fileio->read_file($datafile, \@isolates);
	my $num_lines = scalar @isolates;
	unless ($num_lines) {
		print "\n\t # No data read from file '$datafile'\n\n";
		next;		
	}
	
	# Offer options for top-level stratification 
	my $uc_species = uc $species;
    my $list_question = "\n\n\t Please select data fields for stratifying $uc_species isolates:";
  	my $choice;
  	my %options;
  	$options{1} = 'Subtype';
  	$options{2} = 'Subtype+Region';
  	$options{3} = 'Subtype+Host';
  	print "\n\t # Option 1: Subtype";
	print "\n\t # Option 2: Subtype+Region";
	if ($species eq 'iav') {
		print "\n\t # Option 3: Subtype+Host";
		print "\n\t # Option 4: Exit";
	    $choice = $console->ask_list_question($list_question, 4);
	    if ( $choice eq 4 ) { print "\n\t Exiting...\n\n\n"; exit; }
	}
	else {
		print "\n\t # Option 3: Exit";
	    $choice = $console->ask_list_question($list_question, 3);
	    if ( $choice eq 3 ) { print "\n\t Exiting...\n\n\n"; exit; }
	}

    # Set a minimum number of isolates
    # TODO
    my $min = 1; 
    
    # Set a maximum number of isolates
    # TODO
    my $max = undef; 

	# Stratify isolates based on the user choice
	my %stratified;
	stratify_isolates(\@isolates, \%options, $choice, \%stratified, $min, $max);
	#$devtools->print_hash(\%stratified); die;

	# Show how the data will break down
	my @selected;
	select_isolates(\%stratified, \@selected);
	#$devtools->print_array(\@selected); die;

	# create GLUE program logic for downloading complete genome isolates
	foreach my $key (@selected) { 	# Iterate through export sets

		my $isolates_array_ref =  $stratified{$key};
		#$devtools->print_hash($isolates_array_ref); die;
      
		# Key needs to be used in file names so deal with any char-related issues		
		$key =~ s/:/_/g; # Swap colon for hyphen
		$key =~ s/ /_/g; # Remove any white space (e.g. in host species or country name)
		my $lc_key = lc $key;
		$lc_key =~ s/_/-/g;
		my $descriptorG = 'ncbi-nuccore-' . $lc_key; # For source name 
		my $camel_case_key = $key;
		$camel_case_key =~ s/_([a-z])/\U$1/g; # Convert to camel case
		$camel_case_key =~ s/_//g; # Remove underscore
		my $descriptorM = 'NcbiNuccore' . $camel_case_key;  # For module name
		print "\n\t Key:\t$key";
		print "\n\t Source:\t$descriptorG";
		print "\n\t Module:\t$descriptorM";

		# Write isolate-ordered tabular data for complete genome sequences
		
		my $outfile  = 	$descriptorG . '.tsv';
		write_isolate_data_table($isolates_array_ref, $outfile, $species); 
		create_glue_download_program_logic($isolates_array_ref, $species, $descriptorG, $descriptorM);
	
	}

	print "\n\t Exiting...\n\n\n"; exit;

}

#***************************************************************************
# Subroutine:  select_isolates
# Description: select isolate subsets for export and downloading
#***************************************************************************
sub select_isolates {

	my ($stratified_ref, $selected_ref) = @_;

	# Get user to select subsets
	my %selection_sets = %$stratified_ref;	
	# Iterate
	my $done = undef;
	do {

		# Show stratification summary
		my %select;
		my $num_options = show_stratification_summary(\%selection_sets, \%select); 
		$num_options++;
		my $add_all_option = $num_options;
		print "\n\t # $add_all_option\t\tAdd all and exit";
		$num_options++;
		my $exit_option = $num_options;
		print "\n\t # $exit_option\t\tExit";
		
		my $question .= "\n\n\t # Enter a number to select sets, or to exit selection:";
		my $selected = $console->ask_list_question($question, $num_options);
		
		if ($selected eq $add_all_option) {
		
			# Add all and exit;
			my @selection_keys = keys %select;
			foreach my $selection_key (@selection_keys) {
				
				my $set_key = $select{$selection_key};
				push (@$selected_ref, $set_key);
	
			}
			$done = 1;			
		}
		elsif ($selected eq $exit_option) {
		
			$done = 1;
			
		}

		else {
		
			# Record the selection
			my $set_key = $select{$selected};						
			push (@$selected_ref, $set_key);

			# Remove the selected element from the selection hash
			delete $selection_sets{$set_key};
		
		}

	} until ($done);
	

}

#***************************************************************************
# Subroutine:  show_stratification_summary
# Description: 
#***************************************************************************
sub show_stratification_summary {

	my ($stratified_ref, $select_ref) = @_;

	my @keys = keys %$stratified_ref;
	my $i;
	my %counts;
	foreach my $key (@keys) {
				
		$i++;
		my $layer_isolates_ref = $stratified_ref->{$key};
		my $num_isolates = scalar keys %$layer_isolates_ref;
		$counts{$key} = $num_isolates;
	
	}
		
	# Show the numbers of isolates in each category, in a sorted list
	#my @sorted_keys = sort { $counts{$a} <=> $counts{$b} } keys %counts; # Sort by values
	my @sorted_keys;
	sort_hash_keys($stratified_ref, \@sorted_keys); # Sort by stratifier
	my $j;
	print "\n\t # SUMMARY OF STRATIFICATION RESULTS:\n";
	foreach my $key (@sorted_keys) {

		$j++;
		my $num_isolates = $counts{$key};
		print "\n\t # $j\t\t$key: $num_isolates isolates";
		$select_ref->{$j} = $key;
	}

	return $j;
}

#***************************************************************************
# Subroutine:  stratify_isolates
# Description: 
#***************************************************************************
sub stratify_isolates {

	my ($isolates_ref, $options_ref, $choice, $stratified_ref) = @_;

	my %country_region_data;
	my $region_setting;
	if ($options_ref->{$choice} eq 'Subtype+Region') {
		
		get_country_region_data(\%country_region_data);
		#$devtools->print_hash(\%country_region_data); die;		
		$region_setting = get_region_setting(\%country_region_data);
	
	}

	# Set up column headers
	my $header_row = shift @$isolates_ref;
	chomp $header_row;
	my @header_row = split("\t", $header_row);
	my %column_headers;
	my $i=0;
	foreach my $value (@header_row) {	
		$i++;
		$column_headers{$i} = $value;
	}

	# Index the isolate data by stratifier 'keys'
	my $j;
	foreach my $line (@$isolates_ref) {
	
        $j++;
 	    chomp $line;  
 	    my %isolate_data;
		get_row_values_by_column_header($line, \%column_headers, \%isolate_data);
	    #$devtools->print_hash(\%isolate_data); die;

	    my $isolate     = $isolate_data{'isolate'};
	    my $hn_subtype  = $isolate_data{'gb_subtype'};
	    my $iso_country = $isolate_data{'iso_country'};    
	    my $iso_host    = $isolate_data{'iso_host'};

	    if ($hn_subtype eq '-') {
	    	$hn_subtype = 'NK';
	    }	    
	    if ($iso_country eq '-') {
	    	$iso_country = 'NK';
	    }
	    if ($iso_host eq '-') {
	    	$iso_host = 'NK';
	    }
	    
	    my $key = $hn_subtype;
		if ($options_ref->{$choice} eq 'Subtype+Region') {
			my $region = $country_region_data{$iso_country}->{$region_setting};
			unless ($region) { $region = 'NK'; }
			$key .= ':' . $region;
		}
		elsif ($options_ref->{$choice} eq 'Subtype+Host') {
			$key .= ':' . $iso_host;
		}
		
		if ($stratified_ref->{$key}) {
			my $layer_isolates_ref = $stratified_ref->{$key};
			$layer_isolates_ref->{$isolate} = \%isolate_data;		
		}
		else {
			my %layer_isolates;
			$layer_isolates{$isolate} = \%isolate_data;
			$stratified_ref->{$key} = \%layer_isolates;		
		}
	
	}

}

#***************************************************************************
# Subroutine:  get_region_setting
# Description: index the regions by country
#***************************************************************************
sub get_region_setting {

	my ($country_region_data_ref) = @_;

	# Offer options for top-level stratification 
    my $list_question = "\n\n\t Please select granularity for categorisation of regions:";

  	my %options;
  	$options{1} = 'country';
  	$options{2} = 'sub-region';
  	$options{3} = 'region';
  	print "\n\t # Option 1: Country";
	print "\n\t # Option 2: Sub-region";
	print "\n\t # Option 3: Region";
    my $choice = $console->ask_list_question($list_question, 3);  
    my $setting = $options{$choice};
    
    return $setting;  
 		
}

#***************************************************************************
# Subroutine:  get_country_region_data
# Description: index the regions by country
#***************************************************************************
sub get_country_region_data {

	my ($data_ref) = @_;

    # Read in the file with ISO region data
    my @data_file;
    unless ($fileio->read_file($country_region_data_path, \@data_file)) {
        die "\n\t # Error reading file '$country_region_data_path'\n\n";
    }
	
    # Set up column headers
    my $header_row = shift @data_file;
    chomp $header_row;
    my @header_columns = split("\t", $header_row);

	# Index the data by country
	foreach my $line (@data_file) {
	
        #chomp $line;
        my %result;
        
        @result{@header_columns} = split("\t", $line);
		#$devtools->print_hash(\%result); die;
        my $country  = $result{'name'};
        my $region  = $result{'region'};
        my $subregion  = $result{'sub-region'};

		if ($data_ref->{$country}) {

			$data_ref->{'country'} = $country;
			if ($region) { 
				$data_ref->{'region'} = $region;
			}
			if ($subregion) {
				$data_ref->{'sub-region'} = $subregion;
			}
		
		}
		else {
			my %country_data;
			$country_data{'country'} = $country;
			if ($region) {
				$country_data{'region'} = $region; 
			}
			else {
				die;
			}
			if ($subregion) { 
				$country_data{'sub-region'} = $subregion; 
			}
			$data_ref->{$country} = \%country_data;
						
		}
	
	}

}

#***************************************************************************
# Subroutine:  order_by_isolate
# Description: Index influenza GenBank entries by strain name
# Arguments: $datafile (array containing GenBank entry table)
#            $ordered_by_isolate_ref (hash to store GenBank entries indexed by strain)
#***************************************************************************
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

#***************************************************************************
# Subroutine:  check_consistency
# Description: check the consistency of isolate-associated data across influenza segments
#***************************************************************************
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

#***************************************************************************
# Subroutine:  check_completeness
# Description: sort isolates into those that are complete versus incomplete
#***************************************************************************
sub check_completeness {

	my ($ordered_by_isolate_ref, $complete_ref, $incomplete_ref, $species) = @_;
		
	# Iterate through isolates
	my @isolate_ids = keys %$ordered_by_isolate_ref;
	foreach my $isolate_id (@isolate_ids) {
	
		my $isolate_ref = $ordered_by_isolate_ref->{$isolate_id};
		#$devtools->print_hash($isolate_ref); die;
	
	    my %segment_presence;
	    foreach my $segment_num (@segments) {
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

#***************************************************************************
# Subroutine:  compress_isolates
# Description: For each isolate, select the longest sequence entry for each segment
# Arguments: $isolates_ref (hash reference containing consistent isolates)
#            $compressed_ref (hash reference to store compressed isolates)
#***************************************************************************
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

#***************************************************************************
# Subroutine:  write_isolate_data_table
# Description: Write isolate data to a tabular file
# Arguments: $isolates_ref (hash reference containing uncompressed isolates data:
#            (i.e. with redundant segment sequences included)
#            $output_file (output file path)
#            $species (species type to determine number of segments)
#***************************************************************************
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

#***************************************************************************
# Subroutine:  export_gb_entries_from_isolate_hash
# Description: export a tabular file of GenBank entry details from isolate hash
#***************************************************************************
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

############################################################################
# SPECIES-SPECIFIC
############################################################################

#***************************************************************************
# Subroutine:  determine_iav_genome_lineage
# Description: work our the 'complete genome' taxonomy of an influenza A virus
#              isolate based on genotyping of each segment
#***************************************************************************
sub determine_iav_genome_lineage {

	my ($isolate_details_ref, $counts_ref) = @_;

	my %isolate_subtypes;
	my %isolate_hn_subtypes;

	my %inconsistent_seg_subtype_counts;
	my %subtype_counts;
	my %cg_hn_subtype_counts;
	my %incomplete_hn_subtype_counts;

	my @isolate_ids = keys %$isolate_details_ref;	
	foreach my $isolate_id (@isolate_ids) {
		
		# Get all sequences from this strain 
		my $isolate_ref = $isolate_details_ref->{$isolate_id};
		my %segment_subtypes;
		foreach my $segment_num (@segments) {	
			$segment_subtypes{$segment_num} = '-';
		}

		my @segment_nums = keys %$isolate_ref;
		my $i;
		my $seq_count;
		foreach my $segment_num (@segment_nums) {

			if ($segment_num eq 'seq_count') {
				next;
			}

			$i++;
			#print "\n\t SEGMENT '$segment_num':";
			
			my $segment_seqs_ref = $isolate_ref->{$segment_num};		
			my @segment_seq_accessions = keys %$segment_seqs_ref;
			#$devtools->print_array(\@segment_seq_accessions); # DEBUG
			foreach my $accession (@segment_seq_accessions) {
				
				my $isolate_seq_details = $isolate_ref->{$segment_num}->{$accession};
				my $rec_subtype =  $isolate_seq_details->{'rec_subtype'};
				
				my $current_subtype = $segment_subtypes{$segment_num};
				if ($current_subtype ne '-') {
				
					# check if consistent
					my $current_subtype = $segment_subtypes{$segment_num};
					
					if ($rec_subtype ne $current_subtype) {
						
						print "\n\t SUBTYPE INCONSISTENT for segment $segment_num of strain $isolate_id";
						print "\n\t\t '$current_subtype' versus '$rec_subtype'";
												
					}		
				}
				else {
				
					$segment_subtypes{$segment_num} = $rec_subtype;
				
				}
				
			}
			
		}
		
		# Define genome subtypes
		my @segments = sort by_number keys %segment_subtypes;
		my @subtypes;
		my $missing_segment = undef;
		my @hn_subtype;
		foreach my $segment_num (@segments) {
		
			my $subtype = $segment_subtypes{$segment_num};
			if ($subtype eq '-') {
				$missing_segment  = 'true';
			}
			
			if ($segment_num eq 4) {
				push(@hn_subtype, $subtype)
			}
			elsif ($segment_num eq 6) {
				push(@hn_subtype, $subtype)
			}
			
			push(@subtypes, $subtype);
		
		}
		my $hn_subtype = join('', @hn_subtype);
		my $cg_subtype = join('|', @subtypes);
		
		$isolate_subtypes{$isolate_id} = $cg_subtype;
		$isolate_hn_subtypes{$isolate_id} = $hn_subtype;

		unless ($missing_segment) {
			if ($subtype_counts{$cg_subtype}) {
				$subtype_counts{$cg_subtype}++;	
			}
			else {
				$subtype_counts{$cg_subtype} = 1;
			}

			if ($cg_hn_subtype_counts{$hn_subtype}) {
				$cg_hn_subtype_counts{$hn_subtype}++;	
			}
			else {
				$cg_hn_subtype_counts{$hn_subtype} = 1;
			}

		}
		else {
		
		
			if ($incomplete_hn_subtype_counts{$hn_subtype}) {
				$incomplete_hn_subtype_counts{$hn_subtype}++;	
			}
			else {
				$incomplete_hn_subtype_counts{$hn_subtype} = 1;
			}
		
		}
	}

	#$devtools->print_hash(\%subtype_counts); #die;
	my $file = 'cg_subtype_counts.out.tsv';
	my @file;
	my $header_line = 'cg_subtype' . "\t" . 'strain_count' . "\n";
	push (@file, $header_line);
	my @cg_subtypes = keys %subtype_counts;
	foreach my $subtype (@cg_subtypes) {
	
		my $count = $subtype_counts{$subtype};
		my $line = $subtype . "\t" . $count . "\n";
		push (@file, $line);
	}	
	$fileio->write_file($file, \@file);


	# H/N subtypes
	#$devtools->print_hash(\%cg_hn_subtype_counts); die;
	my $hn_file = 'cg_hn_subtype_counts.out.tsv';
	my @hn_file;
	my $hn_header_line = 'cg_hn_subtype' . "\t" . 'strain_count' . "\n";
	push (@hn_file, $hn_header_line);
	my @hn_subtypes = keys %cg_hn_subtype_counts;
	foreach my $hn_subtype (@hn_subtypes) {
	
		my $count = $cg_hn_subtype_counts{$hn_subtype};
		my $hn_line = $hn_subtype . "\t" . $count . "\n";
		push (@hn_file, $hn_line);
	}
	$fileio->write_file($hn_file, \@hn_file);
	
	# Store in hash reference
	$counts_ref->{'cg_subtype_counts'} = \%subtype_counts;
	$counts_ref->{'cg_hn_subtype_counts'} = \%cg_hn_subtype_counts;
	$counts_ref->{'incomplete_hn_subtype_counts'} = \%incomplete_hn_subtype_counts;
	$counts_ref->{'inconsistent_seg_subtype_counts'} = \%inconsistent_seg_subtype_counts;


}

############################################################################
# Utilities
############################################################################

#***************************************************************************
# Subroutine:  remove_sequences_from_source_directory
# Description: 
#***************************************************************************
sub remove_sequences_from_source_directory {

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

#***************************************************************************
# Subroutine:  summarise_influenza_sequence_set
# Description: counts based on each individual segment
#***************************************************************************
sub summarise_influenza_sequence_set {

	my ($datafile) = @_;

	my @data_file;
	$fileio->read_file($datafile, \@data_file);
	my $num_lines = scalar @data_file;
	unless ($num_lines) {
		print "\n\t # No data read from file '$datafile'\n\n";
		next;		
	}
	 
	my $header_row = shift @data_file;
	chomp $header_row;
	my @header_row = split("\t", $header_row);
	my %column_headers;
	my $i=0;
	foreach my $value (@header_row) {
	
		$i++;
		$column_headers{$i} = $value;

	}
	#$devtools->print_hash(\%column_headers);

    # Set up hashes for counting fields
	my %gb_segment;
	my %rec_segment;
    my %countries;
	my %host_species;
	my %years;
	my %rec_subtype;
	my %gb_subtype;
	my $k=0;
	
	foreach my $line (@data_file) {

        $k++;
 	    chomp $line;  
 	    my %result;
		get_row_values_by_column_header($line, \%column_headers, \%result);
	    #$devtools->print_hash(\%result); die;


	    my $iso_host = $result{'iso_host'};
	    
	    my $iso_country_long = $result{'iso_country'};    
	    my @iso_country = split(':', $iso_country_long);
	    my $iso_country = shift @iso_country;
	    
	    my $gb_segment = $result{'gb_segment'};   
	    my $rec_segment = $result{'rec_segment'};
	 
	    my $isolate   = $result{'strain'};	    
	    my $gb_subtype = $result{'gb_subtype'};
	    my $rec_subtype = $result{'rec_subtype'};
	    
	    my $iso_year = $result{'iso_year'};

	    #print "\n\t STRAIN: $isolate";
	    #print "\n\t GB subtype: $gb_subtype";
	    #print "\n\t REC subtype: $rec_subtype";
	    #print "\n\t GB segment: $gb_segment";
	    #print "\n\t ISO year: $iso_year";
	    #print "\n\t ISO host: $iso_host";
	    #print "\n\t ISO country: $iso_country";

	 	
        # Host species	 	
	 	if ($host_species{$iso_host}){
	 		$host_species{$iso_host}++;
	 	}
	 	else {
	 		$host_species{$iso_host} = 1;
	 	}

        # Countries
	 	if ($countries{$iso_country}) {
	 		$countries{$iso_country}++;
	 	}
	 	else {
	 		$countries{$iso_country} = 1;
	 	}

        # Isolation years 	 	 	
	 	if ($years{$iso_year}){
	 		$years{$iso_year}++; 
	 	}
	 	else {
	 		$years{$iso_year} = 1;
	 	}

        # GenBank segment 	
	 	if ($gb_segment{$gb_segment}){
	 		$gb_segment{$gb_segment}++; 	
	 	}
	 	else {
	 		$gb_segment{$gb_segment} = 1;
	 	}	 	
	 	 	
        # REC segment 	 	
	 	if ($rec_segment{$rec_segment}){
	 		$rec_segment{$rec_segment}++; 	
	 	}
	 	else {
	 		$rec_segment{$rec_segment} = 1;
	 	}

	}

	# Show results
    print "\n\t ######## Countries\n";
	$devtools->print_hash(\%countries);
	print "\n\t ######## Host species\n";
	$devtools->print_hash(\%host_species);
	print "\n\t ######## GB segments\n";
	$devtools->print_hash(\%gb_segment);
	print "\n\t ######## Recogniser segments\n";
	$devtools->print_hash(\%rec_segment);
	print "\n\t ######## Sampling years\n";
	$devtools->print_hash(\%years);
	#$devtools->print_hash(\%column_headers);

}

#***************************************************************************
# Subroutine:  convert_iav_strain_names_to_data_fields
# Description: 
#***************************************************************************
sub convert_iav_strain_names_to_data_fields {

	my ($datafile) = @_;

	my @data_file;
	$fileio->read_file($datafile, \@data_file);
	my $num_lines = scalar @data_file;
	unless ($num_lines) {
		print "\n\t # No data read from file '$datafile'\n\n";
		next;		
	}
	 
	my $header_row = shift @data_file;
	chomp $header_row;
	my @header_row = split("\t", $header_row);
	my %column_headers;
	my $i=0;
	foreach my $value (@header_row) {
	
		$i++;
		$column_headers{$i} = $value;

	}
	#$devtools->print_hash(\%column_headers);

	my $k=0;
	my @result;
    my @new_header_fields = qw [ s_genotype s_host s_country s_number s_year ];
    my $new_header_fields = join ("\t", @new_header_fields);
	$header_row .= "\t$new_header_fields\n";
	push(@result, $header_row);

	foreach my $line (@data_file) {

        $k++;
 	    chomp $line;  
 	    my %result;
		get_row_values_by_column_header($line, \%column_headers, \%result);
	    #$devtools->print_hash(\%result); die;
	 
	    my $isolate = $result{'strain'};

        #$isolate =~ tr{/}{_};
	
	    #$isolate ~= s/\//_/g;
	    print "\n\t STRAIN: $isolate";

	
	    my @strain = split('/', $isolate);
	    #$devtools->print_array(\@strain); die;

	 	my $n_genotype = 'Undefined'; 	
	 	my $n_host = 'NK';
		my $n_country = 'Undefined';
		my $n_number = 'Undefined';
		my $n_year = 'Undefined';
	    
	 	if (scalar @strain == 5) {

			$n_genotype = shift @strain;	 	
			$n_host = shift @strain;
			$n_country = shift @strain;
			$n_number = shift @strain;
			$n_year = shift @strain;	


		}
	 	elsif (scalar @strain == 4) {
	 	
	 	    $n_genotype = shift @strain;
	 		$n_country = shift @strain;
	 		$n_number = shift @strain;
	 		$n_year = shift @strain;
	 		

	 	}
	 	elsif (scalar @strain > 5 or scalar @strain < 4) {

	 		$n_country = 'not-standard';
	 		$n_number = 'not-standard';
	 		$n_year = 'not-standard';
	 	
			
	 	}
	 	
	 	my @new_fields = ( $n_genotype , $n_host, $n_country, $n_number, $n_year);
	 	#$devtools->print_array(\@new_fields); die;
	 	my $new_fields = join ("\t", @new_fields);
	 	$line .= "\t$new_fields\n";
	 	push(@result, $line);
	 	
	 
	}
	
	my $outfile = $datafile . '.extended.tsv';
	$fileio->write_file($outfile, \@result);	

}

############################################################################
# Compile GLUE program logic for influenza virus isolate importation
############################################################################

#***************************************************************************
# Subroutine:  create_glue_download_program_logic
# Description: 
# Create GLUE download module configs, one per segment, for each isolate in $isolates_ref
# Create the GLUE code file that runs the module function for each segement
#***************************************************************************
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

#***************************************************************************
# Subroutine:  adjust_segments
# Description: adjust segment number based on the influenzavirus species
#***************************************************************************
sub adjust_segments {

	my ($species) = @_;

	# Check species defined
	unless ($species) { die $USAGE; }
	$species = lc $species; # Make lowercase
	
	unless ($species eq 'iav' or $species eq 'ibv' or $species eq 'icv' or $species eq 'idv') {
	    die "\n\t Invalid value for --species. Must be one of: iav, ibv, icv, idv\n\n\n"
	}
	if ($species eq 'icv' or $species eq 'idv') {
		# Remove the last segment from the global array 'segments'
		pop @segments; 
	}
}

############################################################################
# BASE
############################################################################

#***************************************************************************
# Subroutine:  get_row_values_by_column_header
# Description: 
#***************************************************************************
sub get_row_values_by_column_header {

	my ($line, $column_header_ref, $result_ref) = @_;

	my @line = split("\t", $line);
	my $j=0;
	foreach my $value (@line) {

		$j++;
		my $column_name = $column_header_ref->{$j};
		$result_ref->{"$column_name"} = $value;
	}
}

#***************************************************************************
# Subroutine:  sort_hash_keys
# Description: Custom sort function to handle keys like 'H1N1', 'H1N2', etc.
#***************************************************************************
sub sort_hash_keys {
    my ($hash_ref, $sorted_keys_ref) = @_;

    @$sorted_keys_ref = sort {
        # Extract parts of the keys including optional extended data
        my ($a_h, $a_n, $a_ext) = $a =~ /^H(\d*)N(\d*)(?::(.*))?$/;
        my ($b_h, $b_n, $b_ext) = $b =~ /^H(\d*)N(\d*)(?::(.*))?$/;

        # Default undefined numerical parts to 0 for comparison
        $a_h = 0 unless defined $a_h && $a_h ne '';
        $a_n = 0 unless defined $a_n && $a_n ne '';
        $b_h = 0 unless defined $b_h && $b_h ne '';
        $b_n = 0 unless defined $b_n && $b_n ne '';

        # Compare first by H and N numerically, then by extended data lexically
        $a_h <=> $b_h || $a_n <=> $b_n || ($a_ext // '') cmp ($b_ext // '');
    } keys %$hash_ref;
}

#***************************************************************************
# Subroutine:  batch_blast_against_rsl
# Description: Batch BLAST against reference sequence library to find most
#              similar sequences to those in an input FASTA file
#***************************************************************************
sub batch_blast_against_rsl {

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



############################################################################
# Command line title blurb 
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: does what it says 
#***************************************************************************
sub show_title {

	#$console->refresh();
	my $title       = 'Influenza Isolate Filter Tool';
	my $version     = '1.0';
	my $description = 'Filter/Organise GenBank Influenza A Virus (IAV) Entries';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<robjgiff@gmail.com>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

    my $HELP   = "\n\n\t ### IAV Filtering";
		$HELP .= "\n\t ### usage: $0 -m[options] -i[infile] -s=[influenzavirus species]";

		$HELP  .= "\n\n\t ### Main functions\n"; 
		$HELP  .= "\n\t  -m=1  Create comprehensive, isolate-oriented database from GenBank entry file";
		$HELP  .= "\n\t  -m=2  Select stratified subsets of isolates for downloading via GLUE";
		$HELP  .= "\n\t  -m=3  Select a phylogeny comparison set using blast & a reference sequence library";	

		$HELP  .= "\n\n\t ### Utility functions\n"; 
		$HELP  .= "\n\t  -m=4  Summarise gb entry database";
		$HELP  .= "\n\t  -m=5  Convert influenza isolate names to data fields";	
		#$HELP  .= "\n\t  -m=6  Remove selected sequences from a source directory";

	print $HELP;
}

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################