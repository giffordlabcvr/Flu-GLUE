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

# Flu Filter modules
#use FluBlastBatcher;
use FluIO;
use FluIsolateUtils;

############################################################################
# Globals
############################################################################

# Version number
my $program_version = '0.1 beta';

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

#my $batcher  = FluBlastBatcher->new($fileio, $blast_obj, $devtools);

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
	}
	else {
		die $USAGE;
	}
	print "\n\n\t # Exit\n\n\n";
}

############################################################################
# Compile isolate database
############################################################################

###############################################################################
# Function:     create_isolate_database
# Description:  Entry point for building an isolate database from exported
#               GenBank data. Supports isolate grouping, consistency checks,
#               segment completeness, and optional GLUE logic generation.
# Parameters:   $datafile - input TSV file containing GenBank sequence metadata
#               $species  - viral species identifier (e.g. 'iav', 'ibv')
#               $console  - Console object for interactive prompts
###############################################################################
sub create_isolate_database {

    my ($datafile, $species, $console) = @_;

    # Step 1: Adjust the expected segment list based on species (e.g. IAV = 8 segments)
    my $num_segments = get_num_segments($species);
    my $isolate_utils = FluIsolateUtils->new($fileio, $num_segments, $devtools);
    my $fluio  = FluIO->new($fileio,  $num_segments, $module_path, $glue_dl_code_directory, $devtools);

    # Step 2: Organize sequences into isolate groups based on ID
    my (%ordered_by_isolate, %null_isolate_entries);
    $isolate_utils->order_by_isolate($datafile, \%ordered_by_isolate, \%null_isolate_entries);
	#$devtools->print_hash(\%ordered_by_isolate); die;

    # Step 3: Check for metadata consistency within each isolate group
    my (%consistent, %inconsistent);
    $isolate_utils->check_consistency(\%ordered_by_isolate, \%consistent, \%inconsistent);

    # Step 4: Export inconsistent entries to a diagnostic TSV file
    my $outfile_inconsistent = $species . "_inconsistent_isolate_entries.tsv";
    $fluio->export_gb_entries_from_isolate_hash($outfile_inconsistent, \%inconsistent);

    # Step 5: Assess segment completeness for consistent isolate entries
    my (%complete, %incomplete);
    $isolate_utils->check_completeness(\%consistent, \%complete, \%incomplete, $species);

    # Step 6: Compress segment entries so only one sequence per segment per isolate remains
    my (%compressed_complete, %compressed_incomplete);
    $isolate_utils->compress_isolates(\%complete, \%compressed_complete);
    $isolate_utils->compress_isolates(\%incomplete, \%compressed_incomplete);

    # Step 7 (Optional): For IAV only, determine the complete genome lineage
    if ($species eq 'iav') {
        my %counts;
        # determine_iav_genome_lineage(\%compressed_complete, \%counts); # Currently disabled
    }

    # Step 8: Export tabular summary for complete and incomplete isolates
    my $outfile_complete   = $species . "_complete_isolates.tsv";
    my $outfile_incomplete = $species . "_incomplete_isolates.tsv";
    $fluio->write_isolate_data_table(\%compressed_complete, $outfile_complete, $species);
    $fluio->write_isolate_data_table(\%compressed_incomplete, $outfile_incomplete, $species);

    # Step 9: Optionally generate GLUE code for downloading complete isolates
    my $num_complete_isolates = scalar keys %compressed_complete;
    if ($num_complete_isolates <= 1000) {
        #my $ask = "\n\t  Do you want to export GLUE modules and code for downloading these isolates?";
        #my $answer = $console->ask_yes_no_question($ask);
        #if ($answer eq 'y') {
        #    $fluio->create_glue_download_program_logic(\%compressed_complete, $species, "ncbi-curated", "NcbiCurated");
        #}
    }
}

#***************************************************************************
# Subroutine:  get_num_segments
# Description: adjust number of segments based on the influenzavirus species
#***************************************************************************
sub get_num_segments {

	my ($species) = @_;

	# Check species defined
	unless ($species) { die $USAGE; }
	$species = lc $species; # Make lowercase

	my $num_segments = undef;
	
	if ($species eq 'iav' or $species eq 'ibv') {

		$num_segments = 8;
	}
	elsif ($species eq 'icv' or $species eq 'idv') {

		$num_segments = 7;
	}
	else {	
	    die "\n\t Invalid value for --species. Must be one of: iav, ibv, icv, idv\n\n\n"

	}
	return $num_segments;
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

    my $HELP   = "\n\n\t ### Influenza Sequence Data Filtering Tools";

	$HELP .= "\n\t ### usage: $0 -m[options] -i[infile] -s=[influenzavirus species]";

	$HELP  .= "\n\n\t ### Main Functions\n"; 
	$HELP  .= "\n\t  -m=1  Build a standardized, isolate-oriented database from GenBank-format entries";
	$HELP  .= "\n\t         • Parses and filters raw GenBank data";
	$HELP  .= "\n\t         • Groups sequences by isolate";
	$HELP  .= "\n\t         • Adjusts segment definitions and checks genome completeness";
	$HELP  .= "\n\t         • Outputs curated isolate records for downstream analysis";

	$HELP  .= "\n\t  -m=2  Select stratified isolate subsets for downstream analysis or export";
	$HELP  .= "\n\t         • Applies spatiotemporal and genome coverage filters";
	$HELP  .= "\n\t         • Supports exporting GLUE-ready modules and code for downloading selected isolates";

	$HELP  .= "\n\t  -m=3  Identify most-similar reference sequences via batch BLAST";
	$HELP  .= "\n\t         • Compares input sequences against a reference sequence library (RSL)";
	$HELP  .= "\n\t         • Selects one or more representatives per isolate for tree-building";

	print $HELP;
}

############################################################################
# EOF
############################################################################