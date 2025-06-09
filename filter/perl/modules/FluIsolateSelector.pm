package FluIsolateSelector;

###############################################################################
# Module:       FluIsolateSelector.pm
# Description:  Supports stratification of influenza isolates based on metadata
#               fields (e.g., subtype, host, region), and exports GLUE logic and
#               data tables for selected subsets.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(select_isolate_subsets stratify_isolates select_isolates);

# ------------------ Subroutine Definitions ------------------

sub get_region_setting {

	my ($country_region_data_ref) = @_;

	# Offer options for top-level stratification 
    my $list_question = "

	 Please select granularity for categorisation of regions:";

  	my %options;
  	$options{1} = 'country';
  	$options{2} = 'sub-region';
  	$options{3} = 'region';
  	print "
	 # Option 1: Country";
	print "
	 # Option 2: Sub-region";
	print "
	 # Option 3: Region";
    my $choice = $console->ask_list_question($list_question, 3);  
    my $setting = $options{$choice};
    
    return $setting;  
 		
}

sub get_country_region_data {

	my ($data_ref) = @_;

    # Read in the file with ISO region data
    my @data_file;
    unless ($fileio->read_file($country_region_data_path, \@data_file)) {
        die "
	 # Error reading file '$country_region_data_path'

";
    }
	
    # Set up column headers
    my $header_row = shift @data_file;
    chomp $header_row;
    my @header_columns = split("	", $header_row);

	# Index the data by country
	foreach my $line (@data_file) {
	
        #chomp $line;
        my %result;
        
        @result{@header_columns} = split("	", $line);
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
	print "
	 # SUMMARY OF STRATIFICATION RESULTS:
";
	foreach my $key (@sorted_keys) {

		$j++;
		my $num_isolates = $counts{$key};
		print "
	 # $j		$key: $num_isolates isolates";
		$select_ref->{$j} = $key;
	}

	return $j;
}

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

1;  # End of FluIsolateSelector module
