
var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];

var viruses = [ 'IAV', 'IBV', 'ICV', 'IDV' ];

_.each(viruses, function(virus) {

	_.each(segments, function(segment) {

        if (virus == 'ICV' && segment == 8 || virus == 'IDV' && segment == 8) {                
        	// Do nothing because there is no segment 8 in ICV or IDV
        }
        else {


			var sourceName = virus + '-ncbi-refseqs-seg' + segment;

			// list the sequences in source
			var listSeqResult = glue.command(["list", "sequence", "-w", "source.name = '"+sourceName+"'"]);

			// extract from the result a list of sequence IDs.
			var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");
			glue.log("INFO", "Source name", sourceName);

			// for each sequence ID
			_.each(seqIds, function(seqId) {

				glue.inMode("sequence/"+sourceName+"/"+seqId, function() {

					glue.command(["set", "field", "segment", segment]);				
								
				});

			});
			
		}

	});

});

