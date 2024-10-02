    # Index isolate data by accession number
	my %isolate_by_accession;
	my %segment_by_accession;

	index_isolate_data(\%isolate_by_accession, \%segment_by_accession);
	#$devtools->print_hash(\%isolate_by_accession); die;
	#$devtools->print_hash(\%segment_by_accession); die;

#***************************************************************************
# Subroutine:  index_isolate_data
# Description: generate two hashes that link:
#              (i)  accession to isolate
#              (ii) accession to segment
#***************************************************************************
sub index_isolate_data {

	my ($iso_index_ref, $segment_index_ref) = @_;

	#my $iso_path = $self->{isolate_path};
	my $iso_path;

	my @isolate_data;
	$fileio->read_file($iso_path, \@isolate_data);
	#$devtools->print_array(\@isolate_data); die;
	
	# Set up isolate file column headers
	my %isolate_headers;
	my $isolate_header_row = shift @isolate_data;
	chomp $isolate_header_row;
	my @isolate_header_row = split("\t", $isolate_header_row);	
	my $j=0;
	foreach my $value (@isolate_header_row) {
		$j++;
		$isolate_headers{$j} = $value;
	}

	foreach my $isolate_row (@isolate_data) {
	
		my %isolate_row;
		get_row_values_by_column_header($isolate_row, \%isolate_headers, \%isolate_row);
		#$devtools->print_hash(\%isolate_row); die;

		# Get relevant row values
		my $strain_id = $isolate_row{strain_id};
		my $segment1_accession = $isolate_row{segment1_accession};
		my $segment2_accession = $isolate_row{segment2_accession};
		my $segment3_accession = $isolate_row{segment3_accession};
		my $segment4_accession = $isolate_row{segment4_accession};
		my $segment5_accession = $isolate_row{segment5_accession};
		my $segment6_accession = $isolate_row{segment6_accession};
		my $segment7_accession = $isolate_row{segment7_accession};
		my $segment8_accession = $isolate_row{segment8_accession};

        $segment_index_ref->{$segment1_accession} = '1';
        $segment_index_ref->{$segment2_accession} = '2';
        $segment_index_ref->{$segment3_accession} = '3';
        $segment_index_ref->{$segment4_accession} = '4';
        $segment_index_ref->{$segment5_accession} = '5';
        $segment_index_ref->{$segment6_accession} = '6';
        $segment_index_ref->{$segment7_accession} = '7';
        $segment_index_ref->{$segment8_accession} = '8';

        $iso_index_ref->{$segment1_accession} = \%isolate_row;
        $iso_index_ref->{$segment2_accession} = \%isolate_row;;
        $iso_index_ref->{$segment3_accession} = \%isolate_row;;
        $iso_index_ref->{$segment4_accession} = \%isolate_row;;
        $iso_index_ref->{$segment5_accession} = \%isolate_row;;
        $iso_index_ref->{$segment6_accession} = \%isolate_row;;
        $iso_index_ref->{$segment7_accession} = \%isolate_row;;
        $iso_index_ref->{$segment8_accession} = \%isolate_row;;

	}	

}





        my $i = '0';
        my @closest_blast_matches;
        do {
            
            # code block
            my $hit_ref = $results[$i];
            #$devtools->print_hash($hit_ref); die;
            my $accession = $hit_ref->{'scaffold'};
            #print "\n\t Match $accession";
            my $segment = $segment_by_accession{$accession};
            $i++;

            if ($segment) {
            
				my $isolate_ref = $isolate_by_accession{$accession};
				#$devtools->print_hash($isolate_ref); die;
				my $hn_subtype  = $isolate_ref->{hn_subtype};
				my $iso_year    = $isolate_ref->{iso_year};
				my $iso_country = $isolate_ref->{iso_country};
				my $strain_id   = $isolate_ref->{strain_id};          
				print "\n\t Match $i: $accession (segment $segment): isolate ($strain_id), country ($iso_country), year ($iso_year), subtype ($hn_subtype)";
            	my $line = "$header\t$seq_id\t$accession\t$segment\t$strain_id\t$iso_country\t$iso_year\t$hn_subtype\n";
            	push (@closest_blast_matches, $line);
            	
            	# Copy the sequence to the a dedicated source directory
            	my $source_file_path = '/Users/rob/Sources/flu/iav-ncbi-curated-segment-' . $segment . '/' . $accession . '.xml';
            	my $destination_path = 'output/iav_refseq_selection_segment' . $segment;
		        my $command = "cp $source_file_path $destination_path";
		        print "\n\t command: $command\n\n"; 
		        system "$command";
            	
            }
            else {
            
 				print "\n\t Match $i: $accession NOT FOUND!";
 				$not_found{$accession} = 1;
           
            }
            
        } until($i eq 500);
