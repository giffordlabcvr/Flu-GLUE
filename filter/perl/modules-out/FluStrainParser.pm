package FluStrainParser;

###############################################################################
# Module:       FluStrainParser.pm
# Description:  Parses IAV strain names into structured metadata fields
#               and appends them to a tab-delimited metadata file.
# Author:       Rob J Gifford
# Version:      1.0
###############################################################################

use strict;
use warnings;
use Exporter 'import';

# Exported functions
our @EXPORT_OK = qw(convert_iav_strain_names_to_data_fields);

# ------------------ Subroutine Definitions ------------------

sub convert_iav_strain_names_to_data_fields {
 
    my ($datafile, $fileio, $devtools) = @_;

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

1;  # End of FluStrainParser module
