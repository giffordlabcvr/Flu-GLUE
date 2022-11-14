var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];

_.each(segments, function(segment) {

    var sourceName = 'iav-ncbi-refseqs-seg' + segment;
    
	// list the sequences in source
	var listSeqResult = glue.command(["list", "sequence", "-w", "source.name = '"+sourceName+"'"]);

	// extract from the result a list of sequence IDs.
	var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");
	//glue.log("INFO", "Source name", sourceName);

	// for each sequence ID
	_.each(seqIds, function(seqId) {

		// create an object in the custom table which uses the sequence ID as the row ID.
		glue.command(["create", "custom-table-row", "isolate_data", seqId]);
	
		// associate the corresponding sequence with this object.
		glue.inMode("sequence/"+sourceName+"/"+seqId, function() {
		
			glue.command(["set", "link-target", "isolate_data", "custom-table-row/isolate_data/"+seqId]);
			
		});
	
	});

});
