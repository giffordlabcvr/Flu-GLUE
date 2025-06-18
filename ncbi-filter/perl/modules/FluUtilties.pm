############################################################################
# Utilities
############################################################################

#***************************************************************************
# Subroutine:  remove_sequences_from_source_directory
# Description: 
#***************************************************************************
sub remove_sequences_from_source_directory {

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

}

#***************************************************************************
# Subroutine:  summarise_influenza_sequence_set
# Description: counts based on each individual segment
#***************************************************************************
sub summarise_influenza_sequence_set {

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


#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	
