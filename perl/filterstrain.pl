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
my @segments = qw [ 1 2 3 4 5 6 7 8 ];

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
	my $file2        = undef;
	
	show_title();

	# Read in options using GetOpt::Long
	GetOptions ('help!'               => \$help,
                'version!'            => \$version,
				'mode|m=i'            => \$mode,
				'infile|i=s'          => \$file1,
				'infile2|f=s'          => \$file2,
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
			filter_out_refseqs_from_strains_file($file1, $file2);
		}
		elsif (($mode eq 3) and $file1) {
			count_across_iav_fields($file1);
		}
		elsif (($mode eq 4) and $file1) {
			convert_iav_strain_names_to_data_fields($file1);
		}
		elsif (($mode eq 5) and $file1) {
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
sub filter_out_refseqs_from_strains_file {

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

	die;

}

#***************************************************************************
# Subroutine:  filter_by_strain
# Description: 
#***************************************************************************
sub filter_by_strain {

	my ($datafile) = @_;

	my %strain_details;
	order_by_strain($datafile, \%strain_details);
	#$devtools->print_hash(\%strain_details); die; # DEBUG

	my %strain_subtypes;
	my %strain_hn_subtypes;
	my %inconsistent_subtype_strain_ids;
	my %counts;
	determine_genome_subtypes(\%strain_details, \%strain_subtypes, \%strain_hn_subtypes, \%inconsistent_subtype_strain_ids, \%counts);
	#$devtools->print_hash(\%strain_subtypes);
	
    # Check for consistency & fill missing values
	my %strains; # flatten strain/isolate associated level to a single level
	my %inconsistent;

    check_consistency_across_strain_sequences(\%strain_details, \%strains, \%inconsistent);
    
    # WRITE strain details
    my @strain_ids = keys %strains;
    my @complete_strain_ids;
    my @complete_strains_file;
    my @incomplete_strains_file;
    my @header_fields = qw [ complete seq_count cg_subtype hn_subtype iso_year iso_month
                             iso_day iso_host iso_source iso_country iso_place_name
                             segment1_accession segment2_accession segment3_accession segment4_accession
                             segment5_accession segment6_accession segment7_accession segment8_accession ];

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
    		print "\n\t  Segment $segment accession: $accession:";# die;
    		$strain_ref->{$key} = $accession;
    		
    	}
    	    	
    	# Set the genome-level subtype
    	my $cg_subtype = $strain_subtypes{$strain_id};
    	my $hn_subtype = $strain_hn_subtypes{$strain_id};
    	
    	
    	unless ($hn_subtype eq 'H11N3') { next; }
    	
    	unless ($cg_subtype and $hn_subtype) {
    		die;
    	} 
		$strain_ref->{'cg_subtype'} = $cg_subtype;
		$strain_ref->{'hn_subtype'} = $hn_subtype;
     	
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
	die;

    # Create download logic for complete genome strains    
    my $subtype_counts_ref = $counts{'cg_subtype_counts'};
	my $cg_hn_subtype_counts_ref = $counts{'cg_hn_subtype_counts'};
   
 	my @hn_subtypes = keys %$cg_hn_subtype_counts_ref;
 	#$devtools->print_hash(\%counts); die;

	# GLUE MASTER
	my @import_source_master;
	my @download_source_master;	
	foreach my $hn_subtype (@hn_subtypes) {
	
		my $count = $cg_hn_subtype_counts_ref->{$hn_subtype};
		create_strain_level_download($hn_subtype, \%strain_details, \%strains, \%strain_hn_subtypes);
		my $comment = "\n\n\t# IAV subtype '$hn_subtype' ($count isolates)";	
		push (@import_source_master, $comment);	
		push (@download_source_master, $comment);
		foreach my $segment_num (@segments) {

			my $source_name = "iav-ncbi-subtype-" . $hn_subtype . "-segment-" . $segment_num;
			my $import_line = "\n\t#import source sources/$source_name";
		
		    my $glue_file_path = 'glue/build/genus/iav/download/importNcbiSequencesIav' . $hn_subtype . '.glue';
			my $download_line = "\n\t#run file $glue_file_path ";
			push (@import_source_master, $import_line);	
			push (@download_source_master, $download_line);
			
		}

	}
	
	# EXPORT
	my $import_master = 'importIavSources.glue';
	my $download_master = 'downloadIavSources.glue';
	
	$fileio->write_file($import_master, \@import_source_master);
	$fileio->write_file($download_master, \@download_source_master);

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
# Subroutine:  create_strain_level_download
# Description: 
#***************************************************************************
sub create_strain_level_download {

	my ($subtype, $strain_details_ref, $strain_subtypes_ref, $strain_hn_subtypes_ref) = @_;
	
	# Create a separate download module for each segment, and...
	# Create the GLUE code file that runs the download function for each segement
	#$devtools->print_hash($strain_subtypes_ref); die;
	
	# HEADER & TAIL
	my @delete_text;
	my @create_text;	
	my @glue_text;	

	my %strain_count;
    my @strain_ids = keys %$strain_subtypes_ref;
	foreach my $segment_num (@segments) {

        my $source_name = "iav-ncbi-curated-segment-" . $segment_num;
		my $header = "<!-- NCBI import module for IAV subtype '$subtype' -->\n<ncbiImporter>\n\t<sequenceFormat>GENBANK_XML</sequenceFormat>\n\t<sourceName>$source_name</sourceName>\n\t<specificPrimaryAccessions>";		
		my $tail = "\n\t</specificPrimaryAccessions>\n\t<sequenceIdField>PRIMARY_ACCESSION</sequenceIdField>\n</ncbiImporter>";
		print "\n\t CREATING DOWNLOAD MODULES FOR '$subtype' SEGMENT '$segment_num':";
		my @module_text;
		push (@module_text, $header);

		foreach my $strain_id (@strain_ids) {
	
			# Get all sequences from this strain
			my $strain_ref = $strain_subtypes_ref->{$strain_id};
            my $details_ref = $strain_details_ref->{$strain_id};
			
			unless ($strain_ref->{'complete'} eq 'TRUE') { 
				next;
			}

			if ($strain_ref->{'hn_subtype'}) {
			
				my $hn_subtype = $strain_ref->{'hn_subtype'};
				unless ($hn_subtype eq $subtype) { 
					next;
				}
			}
			else {
			    next;
				#$devtools->print_hash($details_ref); # die;
				#$devtools->print_hash($strain_ref); # die;			
			}
	
			$strain_count{$strain_id} = 1;
			
			#$devtools->print_hash($details_ref); # DEBUG
			#$devtools->print_hash($strain_ref); # DEBUG

			my $segment_seqs_ref = $details_ref->{$segment_num};		
			my @segment_seq_accessions = keys %$segment_seqs_ref;
			#$devtools->print_array(\@segment_seq_accessions); die; # DEBUG
		
			my $line;
			my $num_segment_seqs = scalar @segment_seq_accessions;
			if ($num_segment_seqs eq 1) {
				my $accession = shift @segment_seq_accessions;
				$line = "\n\t\t<primaryAccession>$accession</primaryAccession> <!-- # $strain_id -->";
			
			}
			else {
			
				# TODO - CHOOSE BEST ONE (Longest?) 
				my $accession = shift @segment_seq_accessions;
				$line = "\n\t\t<primaryAccession>$accession</primaryAccession> <!-- # $strain_id -->";
				
			}
			push (@module_text, $line);
			
		}		
		
		push (@module_text, $tail);
		my $module_name = 'ncbiImporter' . $subtype . 'Segment' . $segment_num;
		my $file_name = 'output/modules/' . $module_name . '.xml';
		$fileio->write_file($file_name, \@module_text);

        my $delete_line = "\n\tdelete module $module_name"; 
		push (@delete_text, $delete_line);
        my $create_line = "\n\tcreate module -f modules/build/genus/iav/downloads/$module_name" . ".xml"; 
		push (@create_text, $create_line);

		my $glue_line = "\n\tmodule $module_name import";
		push (@glue_text, $glue_line);

	}

    # ASSEMBLE GLUE FILE and write as output
    my $strain_count = scalar keys %strain_count;
    my $header_comment = "\t# This GLUE script downloads $strain_count complete genome isolates of IAV subtype '$subtype' from NCBI  GenBank\n\n";

    my $delete_comment = "\t# Delete previous modules if they exist";
    my $create_comment = "\t# Create the download modules";

	my @set_up_text = ( "$delete_comment@delete_text\n\n", "$create_comment@create_text\n\n" );	
	unshift (@set_up_text, $header_comment);

    my $glue_comment .= "\t# Download segment sequences - each is downloaded to a separate source";

	my @glue_output = ( "@set_up_text", "$glue_comment@glue_text\n\n" );
	my $glue_file_name = 'output/glue/downloadNcbiSequencesIav' . $subtype . '.glue';
	$fileio->write_file($glue_file_name, \@glue_output);

}

#***************************************************************************
# Subroutine:  determine_genome_subtypes
# Description: 
#***************************************************************************
sub determine_genome_subtypes {

	my ($strain_details_ref, $strain_subtypes_ref, $strain_hn_subtypes_ref, $inconsistent_ref, $counts_ref) = @_;

	my %inconsistent_seg_subtype_counts;
	my %subtype_counts;
	my %cg_hn_subtype_counts;
	my %incomplete_hn_subtype_counts;

	my @strain_ids = keys %$strain_details_ref;	
	foreach my $strain_id (@strain_ids) {
		
		# Get all sequences from this strain 
		my $strain_ref = $strain_details_ref->{$strain_id};
		my %segment_subtypes;
		foreach my $segment_num (@segments) {	
			$segment_subtypes{$segment_num} = '-';
		}

		#$devtools->print_array(\@segment_nums); # DEBUG
		my @segment_nums = keys %$strain_ref;
		my $i;
		my $seq_count;
		foreach my $segment_num (@segment_nums) {

			if ($segment_num eq 'seq_count') {
				next;
			}

			$i++;
			#print "\n\t SEGMENT '$segment_num':";
			
			my $segment_seqs_ref = $strain_ref->{$segment_num};		
			my @segment_seq_accessions = keys %$segment_seqs_ref;
			#$devtools->print_array(\@segment_seq_accessions); # DEBUG
			foreach my $accession (@segment_seq_accessions) {
				
				my $strain_seq_details = $strain_ref->{$segment_num}->{$accession};
				my $rec_subtype =  $strain_seq_details->{'rec_subtype'};
				
				my $current_subtype = $segment_subtypes{$segment_num};
				if ($current_subtype ne '-') {
				
					# check if consistent
					my $current_subtype = $segment_subtypes{$segment_num};
					
					if ($rec_subtype ne $current_subtype) {
						
						print "\n\t SUBTYPE INCONSISTENT for segment $segment_num of strain $strain_id";
						print "\n\t\t '$current_subtype' versus '$rec_subtype'";
						
						# Record it
						$inconsistent_ref->{$accession} = 1;
						
						if ($inconsistent_seg_subtype_counts{$segment_num}) {
						
							$inconsistent_seg_subtype_counts{$segment_num}++;
						}
						else {
							$inconsistent_seg_subtype_counts{$segment_num} = 1;
						
						}						
						
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
		$strain_subtypes_ref->{$strain_id} = $cg_subtype;
		$strain_hn_subtypes_ref->{$strain_id} = $hn_subtype;

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
			    if ($seq_count > 8) {			    
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
		if ($seg_coverage eq 8) {
			$strains_ref->{$strain_id}->{'complete'} = 'TRUE';
		
		}
		else {
			$strains_ref->{$strain_id}->{'complete'} = 'FALSE';
				
		}	
	}
	
	$devtools->print_hash(\%inconsistent_field_counts);
	
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
	
	if ($field eq 'iso_place_name') {
	  if ($value eq 'NULL') {							
		  $value = $recorded_value;
	  }
	  elsif ($recorded_value eq 'NULL') {							
		  $recorded_value = $value;
	  }
	}
	elsif ($field eq 'iso_host') {

	  my @num_parts = split(' ', $value);
	  my $num_parts = scalar @num_parts;

	  my @num_parts_recorded = split(' ', $recorded_value);							
	  my $num_parts_recorded = scalar @num_parts_recorded;

	  if ($num_parts_recorded > 1 and $num_parts > 1) {

		  # Prefer the more specific of the two
		  if ($num_parts_recorded eq 2 and $num_parts eq 2) {

			  # If one contains 'sp.' and the other does not prefer the one without				
			  if ($recorded_value =~ / sp./) {
				  
				  if ($value =~ / sp./) {
					  print "\n\n\t DID NOT RESOLVE INCONSISTENCY in 2 FIELD '$field' ($value versus $recorded_value) in'$strain_id'!!\n";
				      $inconsistent = 1;
				  }
				  else {
				      print "\n\n\t RESOLVED INCONSISTENCY in FIELD '$field' ($value versus $recorded_value) in'$strain_id'!!\n";

					  $strains_ref->{$strain_id}->{$field} = $value;
				  }
				  
			  }					  
		  }
		  elsif ($num_parts eq 3 and $num_parts_recorded < 3) {
              # Probably 'num_parts' has the subspecies info - use this one
		      $strains_ref->{$strain_id}->{$field} = $value;
		  }
	  }	  
	  elsif ($num_parts > 1) {
		  # Prefer more specific - more components usually better
		  # overwrite existing with non-null value
		  $strains_ref->{$strain_id}->{$field} = $value;	
	  }
	  else {
		  # Prefer the longer of the two
		  my $value_length = length $value;
		  my $recorded_value_length = length $recorded_value;

		  if ($value_length > $recorded_value_length) {
			  $strains_ref->{$strain_id}->{$field} = $value;								    
		  }

	  }

	}

	elsif ($field eq 'iso_source') {
	
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
 	
	 	# skip any that didnt get a recogniser segment
	 	if ($rec_segment =~ '-') {
	 		next;
	 	}

	 	# if GB and recogniser segment agree, keep this one	
		if ($gb_segment eq $rec_segment)  {

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
}

#***************************************************************************
# Subroutine:  count_across_iav_fields
# Description: counts based on each individual segment
#***************************************************************************
sub count_across_iav_fields {

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
	 
	    my $strain   = $result{'strain'};	    
	    my $gb_subtype = $result{'gb_subtype'};
	    my $rec_subtype = $result{'rec_subtype'};
	    
	    my $iso_year = $result{'iso_year'};

	    #print "\n\t STRAIN: $strain";
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
	 
	    my $strain = $result{'strain'};

        #$strain =~ tr{/}{_};
	
	    #$strain ~= s/\//_/g;
	    print "\n\t STRAIN: $strain";

	
	    my @strain = split('/', $strain);
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
	my $title       = 'IAV Filtering Tool';
	my $version     = '1.0';
	my $description = 'Filter/Organise GenBank Influenza A Virus (IAV) Entries';
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
	$HELP  .= "\n\t  -m=2   :   count_across_iav_fields";
	$HELP  .= "\n\t  -m=3   :   convert_iav_strain_names_to_data_fields";	
	$HELP  .= "\n\t  -m=4   :   order_by_strain";

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