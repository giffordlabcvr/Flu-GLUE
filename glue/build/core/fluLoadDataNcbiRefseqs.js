
var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];

var species = [ 'iav', 'ibv', 'icv', 'idv' ];

_.each(species, function(virus) {

	_.each(segments, function(segment) {

		var sourceName = virus + '-ncbi-refseqs-seg' + segment;
	
		// list the sequences in source
		var listSeqResult = glue.command(["list", "sequence", "-w", "source.name = '"+sourceName+"'"]);

		// extract from the result a list of sequence IDs.
		var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");
		glue.log("INFO", "Source name", sourceName);
	
		// for each sequence ID
		_.each(seqIds, function(seqId) {

			// create an object in the custom table which uses the sequence ID as the row ID.
			glue.command(["create", "custom-table-row", "isolate", seqId]);
	
			// associate the corresponding sequence with this object.
			glue.inMode("sequence/"+sourceName+"/"+seqId, function() {
		
				glue.command(["set", "link-target", "isolate", "custom-table-row/isolate/"+seqId]);
				glue.log("INFO", "Sequence", seqId);
				
				if (virus == 'iav') {
					glue.command(["set", "field", "name", 'IAV']);				
					glue.command(["set", "field", "genus", 'Alphainfluenzavirus']);
				}
				if (virus == 'ibv') {			
					glue.command(["set", "field", "name", 'IBV']);	
					glue.command(["set", "field", "genus", 'Betainfluenzavirus']);
				}
				if (virus == 'icv') {
					glue.command(["set", "field", "name", 'ICV']);	
					glue.command(["set", "field", "genus", 'Gammainfluenzavirus']);
				}
				if (virus == 'idv') {			
					glue.command(["set", "field", "name", 'IDV']);	
					glue.command(["set", "field", "genus", 'Deltainfluenzavirus']);
				}
				
				
			});

		});

	});

});
