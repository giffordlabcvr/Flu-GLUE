#!/usr/bin/perl -w
############################################################################
# Script:      flu-gb-filter-tool.pl 
# Description: 
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

# Flu Filter modules
use FluIO;
use FluIsolateUtils;

############################################################################
# Globals
############################################################################

# Version number
my $program_version = '0.1 beta';

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
my ($USAGE) = "\n\t  usage: $0 -i=[infile] -s=[species]\n\n\n";

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
	my ($help, $version, $species, $file1, $file2);
	
	# Read in options using GetOpt::Long
	GetOptions(
		'help!'      => \$help,
		'version!'   => \$version,
		'species|s=s'=> \$species,
		'infile|i=s' => \$file1,
	) or die "Error in command line arguments";

	if ($help) { # Show help page
		show_help_page();  
	}
	elsif ($version)  { 
		print "\n\t # flu-gb-filter-tool.pl version $program_version\n\n";  
	}
	elsif ($species) {
	
		show_title();

	    unless ($file1) {	    
	    	print "\n\t  Please supply an input file using the '-i' option\n";
	    	print "\n\t  For more details on how to run, use '-h' (help)\n";
	    	die $USAGE;
	    }
	    
		create_isolate_database($file1, $species);
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

    # Adjust the expected segment list based on species (e.g. IAV = 8 segments)
    my $num_segments = get_num_segments($species);
    my $module_path = '';            # Empty string -> set to null as this is deprecated
    my $glue_dl_code_directory = ''; # Empty string -> set to null as this is deprecated

    # Instantiate our tools 
    my $isolate_utils = FluIsolateUtils->new($fileio, $num_segments, $species, $devtools);
    my $fluio  = FluIO->new($fileio,  $num_segments, $module_path, $glue_dl_code_directory, $devtools);

    # Organize sequences into isolate groups based on ID
    my (%ordered_by_isolate, %null_isolate_entries);
    $isolate_utils->order_by_isolate($datafile, \%ordered_by_isolate, \%null_isolate_entries);
 
    # Check for metadata consistency within each isolate group
    my (%consistent, %inconsistent);
    $isolate_utils->check_consistency(\%ordered_by_isolate, \%consistent, \%inconsistent);

    # Export inconsistent entries to a diagnostic TSV file
    my $outfile_inconsistent = $species . "_inconsistent_isolate_entries.tsv";
    $fluio->export_gb_entries_from_isolate_hash($outfile_inconsistent, \%inconsistent);

    # Assess segment completeness for consistent isolate entries
    my (%complete, %incomplete);
    $isolate_utils->check_completeness(\%consistent, \%complete, \%incomplete, $species);
	#$devtools->print_hash(\%complete); die;
	#$devtools->print_hash(\%incomplete); die;

    # Step 6: Compress segment entries so only one sequence per segment per isolate remains
    my (%compressed_complete, %compressed_incomplete);
    $isolate_utils->compress_isolates(\%complete, \%compressed_complete);
    $isolate_utils->compress_isolates(\%incomplete, \%compressed_incomplete);

    # Step 8: Export tabular summary for complete and incomplete isolates
    my $outfile_complete   = $species . "_complete_isolates.tsv";
    my $outfile_incomplete = $species . "_incomplete_isolates.tsv";
    $fluio->write_isolate_data_table(\%compressed_complete, $outfile_complete, $species);
    $fluio->write_isolate_data_table(\%compressed_incomplete, $outfile_incomplete, $species);

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
	my $description = 'Build isolate databases from Flu GenBank sequence tables';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<robjgiff@gmail.com>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

    my $HELP   = "\n\n\t ### GenBank Influenza Sequence Data Filtering Tool";

    $HELP .= "\n\n\t ### usage: $0 -i=[infile] -s=[species]";
    $HELP .= "\n\n\t   where [species] must be one of: iav, ibv, icv, idv";

    $HELP  .= "\n\n\t ### Builds a standardized, isolate-oriented influenza sequence database from a table of GenBank entries";
    $HELP  .= "\n\t   • Groups sequences by isolate";
    $HELP  .= "\n\t   • Checks consistency within segments of the same isolate";
    $HELP  .= "\n\t   • Checks completeness across all segments of an isolate";
    $HELP  .= "\n\t   • Outputs curated isolate records for downstream analysis";

    print $HELP;
}

############################################################################
# EOF
############################################################################