	// Run the segment recogniser and capture the results in a map
	
	var resultList;
	glue.inMode("module/iavSegmentRecogniser", function( ) {
	
		resultList = glue.tableToObjects(glue.command(["recognise", "sequence",  "-w", "source.name = 'ncbi-iav-nuccore'"]));		
		_.each(resultList, function(resultObj){

			var sequenceID = resultObj.querySequenceId;
			var categoryId = resultObj.categoryId;
			var direction  = resultObj.direction;

			if (categoryId) {
			  		 glue.logInfo("Processing sequence "+sequenceID+" results: "+categoryId+" ("+direction+")");
			
			}		
		});
		
	});

	
	_.each(resultList, function(resultObj){

		var longSequenceID = resultObj.querySequenceId;
		var nameArr = longSequenceID.split('/');
		var sequenceID = nameArr[nameArr.length - 1];
	
		var categoryId = resultObj.categoryId;
		var direction  = resultObj.direction;

		if (categoryId) {	
			  	
			glue.logInfo("  ADDING DATA FOR sequence "+sequenceID+" results: "+categoryId+" ("+direction+")");
			
			
			glue.inMode("sequence/ncbi-iav-nuccore/"+sequenceID, function() {
			    glue.command(["set", "field", "recogniser_segment", categoryId]);
			});

		}		

		else {
					  	
			glue.logInfo("  EXCLUDING sequence "+sequenceID+" results: "+categoryId+" ("+direction+")");

			glue.inMode("sequence/ncbi-iav-nuccore/"+sequenceID, function() {
				glue.command(["set", "field", "excluded", "TRUE"]);
			});

		}		


	});
	
