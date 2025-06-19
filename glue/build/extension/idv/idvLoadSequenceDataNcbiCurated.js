
var segments = [ 1, 2, 3, 4, 5, 6, 7 ];

_.each(segments, function(segment) {

    var sourceName = 'idv-ncbi-curated-segment-' + segment;
    
	// list the sequences in source
	var listSeqResult = glue.command(["list", "sequence", "-w", "source.name = '"+sourceName+"'"]);

	// extract from the result a list of sequence IDs.
	var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");
	//glue.log("INFO", "Source name", sourceName);
	
	// for each sequence ID
	_.each(seqIds, function(seqId) {

		glue.inMode("sequence/"+sourceName+"/"+seqId, function() {
		
			//glue.log("INFO", "Curated Sequence", seqId);
			glue.command(["set", "field", "name", 'IDV']);
			glue.command(["set", "field", "genus", 'Deltainfluenzavirus']);
			
		});

	});

});

