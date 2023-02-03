#!/usr/bin/perl -w
############################################################################
# Script:      iav-filter.pl 
# Description: a collection of tools for working with sequences + data
# History:     Version 1.0 Creation: Rob J Gifford 2014
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

############################################################################
# Globals
############################################################################

# Version number
my $program_version = '1.0 beta';

# Process ID and time - used to create a unique ID for each program run
my $pid  = $$;
my $time = time;

# Create a unique ID for this process
my $process_id   = $pid . '_' . $time;
my $user = $ENV{"USER"};

# STRAIN INFO - properties that should be shared across all strain sequences
my @strain_info = qw [ iso_year iso_month iso_day iso_host iso_source iso_country iso_place_name ];
my @segment_info = qw [ rec_segment gb_segment gb_subtype rec_subtype ];
my @segments = qw [ 1 2 3 4 5 6 7 ];

############################################################################
# Instantiations for program 'classes' (PERL's Object-Oriented Emulation)
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 -[options]\n\n\n";

############################################################################
# Main program
############################################################################

# Run script
$console->refresh();
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

	# Define options
	my $help         = undef;
	my $version      = undef;
	my $mode		 = undef;
	my $file1        = undef;
	
	show_title();

	# Read in options using GetOpt::Long
	GetOptions ('help!'               => \$help,
                'version!'            => \$version,
				'mode|m=i'            => \$mode,
				'infile|i=s'          => \$file1,
	);

	if ($help) { # Show help page
		show_help_page();  
	}
	elsif ($version)  { 
		print "\n\t # GLUE iav-filter.pl version $program_version\n\n";  
	}
	elsif ($mode)  {
	
	    unless ($file1) {
	    
	    	print "\n\t  Please supply an input file for mode '$mode' using the '-i' option\n";
	    	print "\n\t  For more details on how to run, use '-h' (help)\n";
	    	die $USAGE;
	    }

		if (($mode eq 1) and $file1) {
			filter_by_strain($file1);
		}
		elsif (($mode eq 2) and $file1) {
			count_across_iav_fields($file1);
		}
		elsif (($mode eq 3) and $file1) {
			convert_iav_strain_names_to_data_fields($file1);
		}
		elsif (($mode eq 4) and $file1) {
            my %strain_details;
			order_by_strain($file1, \%strain_details);
		}

	}
	else {
		die $USAGE;
	}
	print "\n\n\t # Exit\n\n\n";
}

############################################################################
# Filter
############################################################################

#***************************************************************************
# Subroutine:  filter_by_strain
# Description: 
#***************************************************************************
sub filter_by_strain {

	my ($datafile) = @_;

	my %strain_details;
	order_by_strain($datafile, \%strain_details);
	#$devtools->print_hash(\%strain_details); die; # DEBUG
	
    # Check for consistency & fill missing values
	my %strains; # flatten strain/isolate associated level to a single level
	my %inconsistent;
    check_consistency_across_strain_sequences(\%strain_details, \%strains, \%inconsistent);
    
    # WRITE strain details
    my @strain_ids = keys %strains;
    my @complete_strain_ids;
    my @complete_strains_file;
    my @incomplete_strains_file;
    my @header_fields = qw [ complete seq_count iso_year iso_month
                             iso_day iso_host iso_source iso_country iso_place_name
                             segment1_accession segment2_accession segment3_accession segment4_accession
                             segment5_accession segment6_accession segment7_accession ];

	foreach my $segment_num (@segments) {	
		my $field_name = 'segment_' .  $segment_num . '_counts';
		push (@header_fields, $field_name);
	}
	my @header_row = ( 'strain_id',  @header_fields);
	my $header_row = join("\t", @header_row);
	push (@complete_strains_file, "$header_row\n");
	push (@incomplete_strains_file, "$header_row\n");
	#print "$header_row"; die;
    #$devtools->print_hash(\%strains); die; # DEBUG
    
    foreach my $strain_id (@strain_ids) {
    	
    	if ($inconsistent{$strain_id}) {    	
    		next;
    	}
    	
    	my @data_line_values;
    	push (@data_line_values, $strain_id);
    	
    	my $strain_ref = $strains{$strain_id};
    	my $strain_details_ref = $strain_details{$strain_id};
    	#$devtools->print_hash($strain_details_ref); # DEBUG
    	
    	# Get the segment accessions
    	my %accessions;
    	select_isolate_segments($strain_details_ref, \%accessions);
    	#$devtools->print_hash(\%accessions); die; # DEBUG
    	foreach my $segment (@segments) {
    	
    		my $accession = $accessions{$segment};
    		my $key = 'segment' . $segment . '_accession';
    		if ($key and $accession) {
        		print "\n\t  Segment $segment accession: $accession:";# die;
    			$strain_ref->{$key} = $accession;
    		
    		}
    	}
    	    	
    	# Set the genome-level subtype     	
    	foreach my $field (@header_fields) {
    	
    		my $value = $strain_ref->{$field};
    		if ($value) {
    			push (@data_line_values, $value);
    			
    		}
    		else {
    			push (@data_line_values, '0');
    			
    		}
    		   	
    	}
    	
    	my $data_line = join("\t", @data_line_values);
    	if ($strain_ref->{'complete'} eq 'TRUE') {    		
			push (@complete_strains_file, "$data_line\n");		
		}
		else {
			push (@incomplete_strains_file, "$data_line\n");
		}

    }
    
    # Write strain-level data for strains with complete genomes (8 segments)
    my $cg_strain_info_file = 'strains.complete.out.tsv'; 
    $fileio->write_file($cg_strain_info_file, \@complete_strains_file);

    # Write strain-level data for strains with incomplete genomes (i.e. <8 segments sequenced)
    my $incomplete_strain_info_file = 'strains.incomplete.out.tsv'; 
    $fileio->write_file($incomplete_strain_info_file, \@incomplete_strains_file);
 
}


#***************************************************************************
# Subroutine:  select_isolate_segments
# Description: get longest if >1
#***************************************************************************
sub select_isolate_segments {

	my ($strain_details_ref, $accessions_ref) = @_;

	foreach my $segment (@segments) {
		
		my $segment_details_ref = $strain_details_ref->{$segment};
		#$devtools->print_hash($segment_details); die; # DEBUG
		
		my @accession = keys %$segment_details_ref;
		my $longest_length;
		my $longest_acc;
		foreach my $accession (@accession) {
		
			my $length = $segment_details_ref->{'length'};
			if ($longest_length) {
								
				if ($length > $longest_length) {
					$longest_length = $length;
					$longest_acc = $accession;
				}				
				
			} else {
				$longest_length = $length;
				$longest_acc = $accession;
							
			}
		
		}
		
		$accessions_ref->{$segment} = $longest_acc;
	
	}

}

#***************************************************************************
# Subroutine:  check_consistency_across_strain_sequences
# Description: Iterate through checking for consistency, missing values
#              Fill in missing values
#              Separate out any strains with internal conflicts	
#***************************************************************************
sub check_consistency_across_strain_sequences {

	my ($strain_details_ref, $strains_ref, $inconsistent_ref) = @_;
	
	my %inconsistent_field_counts;	
	my @strain_ids = keys %$strain_details_ref;
	
	foreach my $strain_id (@strain_ids) {

        my %segment_counts;
        my %segment_subtypes;
		foreach my $segment_num (@segments) {	
			$segment_counts{$segment_num} = '0';
		}
		
		# Get all sequences from this strain 
		my $strain_ref = $strain_details_ref->{$strain_id};

		#$devtools->print_array(\@segment_nums); # DEBUG
		my @segment_nums = keys %$strain_ref;
		my $i;
		my $seq_count;
		foreach my $segment_num (@segment_nums) {

			if ($segment_num eq 'seq_count') {
			    $seq_count = $strain_ref->{'seq_count'};
			    if ($seq_count > 7) {			    
			    	print "\n\t STRAIN '$strain_id' has '$seq_count' sequences";
			    }
				next;
			}

			$i++;
			#print "\n\t SEGMENT '$segment_num':";
			
			my $segment_seqs_ref = $strain_ref->{$segment_num};		
			my @segment_seq_accessions = keys %$segment_seqs_ref;
			#$devtools->print_array(\@segment_seq_accessions); # DEBUG
			foreach my $accession (@segment_seq_accessions) {
				
				my $strain_seq_details = $strain_ref->{$segment_num}->{$accession};
				$segment_counts{$segment_num} = $segment_counts{$segment_num} + 1;
				
				# Iterate through each strain data field, checking for consistency
				foreach my $field (@strain_info) {
			
					#print "\n\t Assessing strain '$strain_id' for consistency in field '$field'";
				
					my $value = $strain_seq_details->{$field};
					my $recorded_value = $strains_ref->{$strain_id}->{$field};
					if ($value and $recorded_value) {
				
						if ($recorded_value eq $value) {
							# consistent, do nothing
						}
						else {
						
							if ($value eq '-') {
								# prefer existing value, do nothing					
							}
	
							elsif ($recorded_value eq '-') {
								# overwrite existing with non-null value
								$strains_ref->{$strain_id}->{$field} = $value;						
							}
							else {
	                           
								my $inconsistent = resolve_conflicts($strain_id, $field, $strains_ref, $strain_seq_details);
								
								
								if ($inconsistent) {
									print "\n\n\t INCONSISTENCY in FIELD ($field) for strain '$strain_id'!!";
	 								print "\n\t VALUE IN HASH = '$recorded_value'";
	                                print "VALUE IN SEQ = '$value'\n";
									$inconsistent_ref->{$strain_id} = 'true';

									if ($inconsistent_field_counts{$field}) {
										$inconsistent_field_counts{$field}++;
								
									}
									else {
										$inconsistent_field_counts{$field} = 1;
									}
									
								}
							}						
						}
					}
					
					elsif ($value) {
						# no recorded value yet, so set one from this sequence				
						$strains_ref->{$strain_id}->{$field} = $value;
					}
				}		
			}
		}
		
		# Add sequence count data
		$strains_ref->{$strain_id}->{'seq_count'} = $seq_count;
		
		my %segment_coverage;
		foreach my $segment_num (@segments) {
		
			my $segment_count = $segment_counts{$segment_num};
			my $field_name = 'segment_' .  $segment_num . '_counts';
			
			$strains_ref->{$strain_id}->{$field_name} = $segment_count;
			if ($segment_count) {
				$segment_coverage{$field_name} = $segment_count;
			}
			
		}
		my @seg_coverage = keys %segment_coverage;
		my $seg_coverage = scalar @seg_coverage;
		if ($seg_coverage eq 7) {
			$strains_ref->{$strain_id}->{'complete'} = 'TRUE';
		
		}
		else {
			$strains_ref->{$strain_id}->{'complete'} = 'FALSE';
				
		}	
	}
	
	#$devtools->print_hash(\%inconsistent_field_counts);
	
}

#***************************************************************************
# Subroutine:  resolve_conflicts
# Description: 
#***************************************************************************
sub resolve_conflicts {

	my ($strain_id, $field, $strains_ref, $strain_seq_details)= @_;

	my $value = $strain_seq_details->{$field};
	my $recorded_value = $strains_ref->{$strain_id}->{$field};

	my $inconsistent;
	
	if ($field eq 'iso_source') {
	
		# Prefer the longer of the two
		my $value_length = length $value;
		my $recorded_value_length = length $recorded_value;
		# Prefer the longer of the two - don't mark as inconsistent as too much noise in this field
		if ($value_length > $recorded_value_length) {
			$strains_ref->{$strain_id}->{$field} = $value;								    
		}

	}

	elsif ($field eq 'iso_day') {

		# Get recorded values for 'iso_month' and 'iso_year'
		# If these don't disagree between the isolates, just take lowest
		my $recorded_year = $strains_ref->{$strain_id}->{'iso_year'}; 
		my $recorded_month = $strains_ref->{$strain_id}->{'iso_month'};
		my $year = $strain_seq_details->{'iso_year'};
		my $month = $strain_seq_details->{'iso_month'};
		
		if ($year eq $recorded_year and $month eq $recorded_month) {
		
			#print "\n\n\t RESOLVED INCONSISTENCY in FIELD ($field) for strain '$strain_id'!!\n";
			if ($value < $recorded_value) {			
				$strains_ref->{$strain_id}->{$field} = $value;			
			}
		}
		else {
		
			$inconsistent = 1;
		}

	}


	else {
	
	  $inconsistent = 1;
		  
	}

	return $inconsistent;	
}

#***************************************************************************
# Subroutine:  order_by_strain
# Description: 
#***************************************************************************
sub order_by_strain {

	my ($datafile, $strain_details_ref) = @_;

	my @data_file;
	$fileio->read_file($datafile, \@data_file);
	my $num_lines = scalar @data_file;
	unless ($num_lines) {
		print "\n\t # No data read from file '$datafile'\n\n";
		next;		
	}
	
	# Set up column headers
	my $header_row = shift @data_file;
	chomp $header_row;
	my @header_row = split("\t", $header_row);
	my %column_headers;
	my $i=0;
	foreach my $value (@header_row) {	
		$i++;
		$column_headers{$i} = $value;
	}

    # Set up hashes for counting fields
	my %rec_subtype;
	my %gb_subtype;
	my $k=0;
	
	foreach my $line (@data_file) {

        $k++;
 	    chomp $line;  
 	    my %result;
		get_row_values_by_column_header($line, \%column_headers, \%result);
	    #$devtools->print_hash(\%result); die;

	    my $sequenceID  = $result{'sequenceID'};	    
	    my $strain      = $result{'strain'};	    
	    my $gb_subtype  = $result{'gb_subtype'};
	    my $rec_subtype = $result{'rec_subtype'};
	    my $gb_segment  = $result{'gb_segment'};
	    my $rec_segment = $result{'rec_segment'};
	    my $lab_host    = $result{'lab_host'};
	    my $iso_host    = $result{'iso_host'};
	    my $iso_country = $result{'iso_country'};
	    my $iso_source  = $result{'iso_source'};
	    my $iso_day     = $result{'iso_day'};
	    my $iso_month   = $result{'iso_month'};
	    my $iso_year    = $result{'iso_year'};
	    my $pubmed_id   = $result{'pubmed_id'};
	    
	    # Ignore NULL strain sequences 
	    if ($strain eq '-') { next; }
	    
	    # Deal with any country field values that also contain place names
	    my @country = split(':', $iso_country);
	    my $array_size = scalar @country;
	    if ($array_size > 1) {
	    	my $iso_place_name = pop @country;
	    	my $iso_country_only = shift @country;
	    	$iso_place_name =~ s/^\s+//; # Remove leading whitespace
	    	
	    	#print "\n\t SPLIT '$iso_country' to: '$iso_country_only' and '$iso_place_name'";	    	
	    	$result{'iso_country'} = $iso_country_only;
	    	$result{'iso_place_name'}  = $iso_place_name; 	
	    	
	    }
	    else {
	    	$result{'iso_place_name'}  = 'NULL'; 	
	    }
 	

		#print "\n\t STRAIN: $strain";
		#print "\n\t SEG: $gb_segment";
		#print "\n\t GB subtype: $gb_subtype";
		#print "\n\t REC subtype: $rec_subtype";

		# Strain details	 	
		if ($strain_details_ref->{$strain}){

			my $strain_ref = $strain_details_ref->{$strain};
			#$devtools->print_hash($strain_ref);
			$strain_ref->{'seq_count'}++;
											
																				 
			if ($strain_ref->{$gb_segment}){
								   
				
				my $strain_segment_seq_details_ref = $strain_ref->{$gb_segment};
				$strain_segment_seq_details_ref->{$sequenceID} = \%result;
																																												
			}
			else {

				my %strain_segment_seq_details;
				$strain_segment_seq_details{$sequenceID} = \%result;
				$strain_ref->{$gb_segment} = \%strain_segment_seq_details;
								
			}	 	

		}
		else {

			my %strain;
			$strain{'seq_count'} = 1;
			
			my %strain_segment_seq_details;
			$strain_segment_seq_details{$sequenceID} = \%result;
			
			$strain{$gb_segment} = \%strain_segment_seq_details;				
			#$strain{$gb_segment} = 1;
			$strain_details_ref->{$strain} = \%strain;
		
		}	 	

		##### Subtype comparisons

		my $gb_segment_subtype = $gb_segment . '-' . $gb_subtype;
		my $rec_segment_subtype = $gb_segment . '-' . $rec_subtype;

		# REC Subtype  	
		if ($rec_subtype{$rec_segment_subtype}){
			$rec_subtype{$rec_segment_subtype}++; 	
		}
		else {
			$rec_subtype{$rec_segment_subtype} = 1;
		}	 	

		# GenBank Subtype 	
		if ($gb_subtype{$gb_segment_subtype}){
			$gb_subtype{$gb_segment_subtype}++; 	
		}
		else {
			$gb_subtype{$gb_segment_subtype} = 1;
		}	 	
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

############################################################################
# Command line title blurb 
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: does what it says 
#***************************************************************************
sub show_title {

	#$console->refresh();
	my $title       = 'ICV and IDV Filtering Tool';
	my $version     = '1.0';
	my $description = 'Filter/Organise GenBank files by strain';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

	my ($HELP) = "\n\t usage: $0 -m[options] -i[infile]";

	$HELP  .= "\n\n\t # IAV Filtering\n";
	$HELP  .= "\n\t  -m=1   :   filter_by_strain";

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

############################################################################