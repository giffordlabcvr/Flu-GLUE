#!/usr/bin/perl -w
############################################################################
# Script:      flu-gb-filter-tool.pl
# Description: A modular toolkit for processing, filtering, annotating, and
#              summarizing influenza virus sequence data and metadata,
#              designed to support isolate-level database construction and
#              downstream genomic analysis workflows.
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

# Refactored modules
use FluStrainParser qw(convert_iav_strain_names_to_data_fields);
use FluSummary qw(summarise_influenza_sequence_set);
use FluSequenceRemover qw(remove_sequences_from_source_directory);
use FluBlastBatcher qw(batch_blast_against_rsl);
use FluIsolateSelector qw(select_isolate_subsets stratify_isolates select_isolates get_region_setting get_country_region_data);
use FluIsolateBuilder qw(create_isolate_database);
use SegmentCoverage qw(adjust_segments check_completeness);
use IsolateUtils qw(order_by_isolate check_consistency compress_isolates);
use FluIO qw(write_isolate_data_table export_gb_entries_from_isolate_hash);
use GlueCodeExporter qw(create_glue_download_program_logic);

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
############################################################################
# Entry point: Main program execution
############################################################################

# Run script
main();

# Exit program
exit;

############################################################################
############################################################################
# Subroutine Definitions
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
            else { print "
      Unknown mode \'$mode\'
";
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
# Subroutine create_isolate_database moved to module


#***************************************************************************
# Subroutine:  select_isolate_subsets
# Description: select subsets of isolate sequences, stratified by their
#***************************************************************************
# Subroutine select_isolate_subsets moved to FluIsolateSelector module


#***************************************************************************
# Subroutine:  select_isolates
# Description: select isolate subsets for export and downloading
#***************************************************************************
# Subroutine select_isolates moved to FluIsolateSelector module


#***************************************************************************
# Subroutine:  show_stratification_summary
# Description:
#***************************************************************************
# Subroutine show_stratification_summary moved to FluIsolateSelector module


#***************************************************************************
# Subroutine:  stratify_isolates
# Description:
#***************************************************************************
# Subroutine stratify_isolates moved to FluIsolateSelector module


#***************************************************************************
# Subroutine:  get_region_setting
# Description: index the regions by country
#***************************************************************************
# Subroutine get_region_setting moved to FluIsolateSelector module


#***************************************************************************
# Subroutine:  get_country_region_data
# Description: index the regions by country
#***************************************************************************
# Subroutine get_country_region_data moved to FluIsolateSelector module


#***************************************************************************
# Subroutine:  order_by_isolate
# Description: Index influenza GenBank entries by strain name
# Arguments: $datafile (array containing GenBank entry table)
#            $ordered_by_isolate_ref (hash to store GenBank entries indexed by strain)
#***************************************************************************
# Subroutine order_by_isolate moved to module


#***************************************************************************
# Subroutine:  check_consistency
# Description: check the consistency of isolate-associated data across influenza segments
#***************************************************************************
# Subroutine check_consistency moved to module


#***************************************************************************
# Subroutine:  check_completeness
# Description: sort isolates into those that are complete versus incomplete
#***************************************************************************
# Subroutine check_completeness moved to module


#***************************************************************************
# Subroutine:  compress_isolates
# Description: For each isolate, select the longest sequence entry for each segment
# Arguments: $isolates_ref (hash reference containing consistent isolates)
#            $compressed_ref (hash reference to store compressed isolates)
#***************************************************************************
# Subroutine compress_isolates moved to module


#***************************************************************************
# Subroutine:  write_isolate_data_table
# Description: Write isolate data to a tabular file
# Arguments: $isolates_ref (hash reference containing uncompressed isolates data:
#            (i.e. with redundant segment sequences included)
#            $output_file (output file path)
#            $species (species type to determine number of segments)
#***************************************************************************
# Subroutine write_isolate_data_table moved to module


#***************************************************************************
# Subroutine:  export_gb_entries_from_isolate_hash
# Description: export a tabular file of GenBank entry details from isolate hash
#***************************************************************************
# Subroutine export_gb_entries_from_isolate_hash moved to module


############################################################################
# SPECIES-SPECIFIC
############################################################################

#***************************************************************************
# Subroutine:  determine_iav_genome_lineage
# Description: work our the 'complete genome' taxonomy of an influenza A virus
#              isolate based on genotyping of each segment
#***************************************************************************
# Subroutine determine_iav_genome_lineage removed (deprecated)


############################################################################
# Utilities
############################################################################

#***************************************************************************
# Subroutine:  remove_sequences_from_source_directory
# Description:
#***************************************************************************
# Subroutine remove_sequences_from_source_directory moved to FluSequenceRemover module


#***************************************************************************
# Subroutine:  summarise_influenza_sequence_set
# Description: counts based on each individual segment
#***************************************************************************
# Subroutine summarise_influenza_sequence_set moved to FluSummary module


#***************************************************************************
# Subroutine:  convert_iav_strain_names_to_data_fields
# Description:
#***************************************************************************
# Subroutine convert_iav_strain_names_to_data_fields moved to FluStrainParser module


############################################################################
# Compile GLUE program logic for influenza virus isolate importation
############################################################################

#***************************************************************************
# Subroutine:  create_glue_download_program_logic
# Description:
# Create GLUE download module configs, one per segment, for each isolate in $isolates_ref
# Create the GLUE code file that runs the module function for each segement
#***************************************************************************
# Subroutine create_glue_download_program_logic moved to module


#***************************************************************************
# Subroutine:  adjust_segments
# Description: adjust segment number based on the influenzavirus species
#***************************************************************************
# Subroutine adjust_segments moved to module


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
# Subroutine batch_blast_against_rsl moved to FluBlastBatcher module




############################################################################
# Command line title blurb
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: does what it says
#***************************************************************************
sub show_title {
    my $title       = 'Influenza Isolate Filter Tool';
    my $version     = '1.0';
    my $description = 'Filter/Organise GenBank Influenza A Virus (IAV) Entries';
    my $author      = 'Robert J. Gifford';
    my $contact        = '<robjgiff@gmail.com>';
    $console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

    my $HELP   = "\n\n\t ### Influenza Virus Filtering";
        $HELP .= "\n\t ### usage: $0 -m[options] -i[infile] -s=[influenzavirus species]";

        $HELP  .= "\n\n\t ### Main functions\n";
        $HELP  .= "\n\t  -m=1  Create comprehensive, isolate-oriented database from GenBank entry file";
        $HELP  .= "\n\t  -m=2  Select stratified subsets of isolates for downloading via GLUE";
        $HELP  .= "\n\t  -m=3  Select a phylogeny comparison set using blast and a reference sequence library";

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
