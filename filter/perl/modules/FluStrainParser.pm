#!/usr/bin/perl -w
############################################################################
# Module:      FluStrainParser.pm 
# Description: 
# History:     Rob Gifford January 2024: Creation
############################################################################
package FluStrainParser;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

############################################################################
# Globals
############################################################################

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new FluStrainParser.pm 'object'
#***************************************************************************
sub new {

    my ($class, $fileio, $devtools) = @_;

	# Member variables
    my $self = {
        fileio    => $fileio,
        devtools  => $devtools,
    };
    bless $self, $class;
    return $self;
}

############################################################################
# Subroutine Definitions
############################################################################

#***************************************************************************
# Subroutine:  convert_iav_strain_names_to_data_fields
# Description: 
#***************************************************************************
sub convert_iav_strain_names_to_data_fields {

	my ($self, $datafile) = @_;

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
1;  # End of FluStrainParser module

