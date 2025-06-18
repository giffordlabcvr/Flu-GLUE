
var seqObjs = glue.tableToObjects(
glue.command(["list", "sequence", "-s", "cluster", "-w", "source.name = 'ncbi-refseqs'", 
			 "sequenceID", "gb_length", "major_clade", "minor_clade", "cluster", "collection_year", "gb_create_date"]));

