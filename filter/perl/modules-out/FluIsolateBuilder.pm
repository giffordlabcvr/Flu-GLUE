package FluIsolateBuilder;

###############################################################################
# Module:       FluIsolateBuilder.pm
# Description:  Builds isolate-oriented datasets from raw GenBank export files.
#               This module is part of the Flu-GLUE filtering pipeline. It groups
#               sequences by isolate ID, checks for metadata consistency and
#               segment completeness, and exports summary tables. Optionally,
#               it generates GLUE download logic for complete isolates.
# Author:       Rob J Gifford
# Version:      1.0
# License:      Open source (license details to be provided in project README)
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Export the main function
our @EXPORT_OK = qw(create_isolate_database);

# Dependencies assumed to be available in the calling script or injected.
# use DevTools;
# use FileIO;
# use Console;

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
    adjust_segments($species);

    # Step 2: Organize sequences into isolate groups based on ID
    my (%ordered_by_isolate, %null_isolate_entries);
    order_by_isolate($datafile, \%ordered_by_isolate, \%null_isolate_entries);

    # Step 3: Check for metadata consistency within each isolate group
    my (%consistent, %inconsistent);
    check_consistency(\%ordered_by_isolate, \%consistent, \%inconsistent);

    # Step 4: Export inconsistent entries to a diagnostic TSV file
    my $outfile_inconsistent = $species . "_inconsistent_isolate_entries.tsv";
    export_gb_entries_from_isolate_hash($outfile_inconsistent, \%inconsistent);

    # Step 5: Assess segment completeness for consistent isolate entries
    my (%complete, %incomplete);
    check_completeness(\%consistent, \%complete, \%incomplete, $species);

    # Step 6: Compress segment entries so only one sequence per segment per isolate remains
    my (%compressed_complete, %compressed_incomplete);
    compress_isolates(\%complete, \%compressed_complete);
    compress_isolates(\%incomplete, \%compressed_incomplete);

    # Step 7 (Optional): For IAV only, determine the complete genome lineage
    if ($species eq 'iav') {
        my %counts;
        # determine_iav_genome_lineage(\%compressed_complete, \%counts); # Currently disabled
    }

    # Step 8: Export tabular summary for complete and incomplete isolates
    my $outfile_complete   = $species . "_complete_isolates.tsv";
    my $outfile_incomplete = $species . "_incomplete_isolates.tsv";
    write_isolate_data_table(\%compressed_complete, $outfile_complete, $species);
    write_isolate_data_table(\%compressed_incomplete, $outfile_incomplete, $species);

    # Step 9: Optionally generate GLUE code for downloading complete isolates
    my $num_complete_isolates = scalar keys %compressed_complete;
    if ($num_complete_isolates <= 1000) {
        my $ask = "\n\t  Do you want to export GLUE modules and code for downloading these isolates?";
        my $answer = $console->ask_yes_no_question($ask);
        if ($answer eq 'y') {
            create_glue_download_program_logic(\%compressed_complete, $species, "ncbi-curated", "NcbiCurated");
        }
    }
}

1; # End of module
