package FluSummary;

###############################################################################
# Module:       FluSummary.pm
# Description:  Computes and prints summary statistics from an influenza
#               sequence metadata table, including host, region, year, segment.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(summarise_influenza_sequence_set);

# ------------------ Subroutine Definitions ------------------

sub summarise_influenza_sequence_set {
    my ($datafile, $fileio, $devtools) = @_;

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

1;  # End of FluSummary module
