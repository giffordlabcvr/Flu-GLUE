// list the sequences in source ncbi-refseqs
var listSourceResult = glue.command(["list", "source"]);

// extract from the result a list of sources
var sourceNames = glue.getTableColumn(listSourceResult, "name");

// for each source
_.each(sourceNames, function(sourceName) {

	//glue.log("INFO", "Source", sourceName);

    var whereClause = "source.name = '"+sourceName+"'";
    
	var listSeqResult = glue.command(["list", "sequence", "-w", whereClause]);

	// extract from the result a list of sequence IDs.
	var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");

	// for each sequence ID
	_.each(seqIds, function(seqId) {

		//glue.log("INFO", "Seqeunce", seqId);

		// create an object in the custom table which uses the sequence ID as the row ID.
		glue.command(["create", "custom-table-row", "isolate", seqId]);
	
		// associate the corresponding sequence with this object.
		glue.inMode("sequence/"+sourceName+"/"+seqId, function() {
			glue.command(["set", "link-target", "isolate", "custom-table-row/isolate/"+seqId]);

		});
		
	});

});
