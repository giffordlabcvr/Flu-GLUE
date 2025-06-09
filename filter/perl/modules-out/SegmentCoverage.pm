package SegmentCoverage;

###############################################################################
# Module:       SegmentCoverage.pm
# Description:  Contains logic for determining the expected segment list based
#               on influenza virus species and evaluating genome completeness
#               for isolate records.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(adjust_segments check_completeness);

# ------------------ Subroutine Definitions ------------------

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

1;  # End of SegmentCoverage module
