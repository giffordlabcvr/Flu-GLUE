
var ncbiNuccoreSeg4;
var whereClause = "source.name like 'ibv-ncbi-nuccore-%' and rec_segment = 4";
ncbiNuccoreSeg4 = glue.tableToObjects(glue.command(["list", "sequence", "source.name", "sequenceID", "-w", whereClause]));
//glue.log("INFO", "RESULT WAS ", ncbiNuccoreSeg4);

var processed = 0;

_.each(ncbiNuccoreSeg4, function(ncbiNuccoreSeg4) {

	var sequenceID = ncbiNuccoreSeg4.sequenceID;
	//glue.log("INFO", "ID RESULT WAS ", sequenceID);

	var sourceName = ncbiNuccoreSeg4["source.name"];
	//glue.log("INFO", "SOURCE RESULT WAS ", sourceName);

	var whereClause = "sequenceID = '" + sequenceID + "'";
	var subtypeResults;
	glue.inMode("/module/ibvSeg4MaxLikelihoodGenotyper", function() {
		subtypeResults = glue.command(["genotype", "sequence", "-w", whereClause]);
		glue.log("INFO", "ibv subtype RESULT WAS ", subtypeResults);			
	});

	var subtypeRows = subtypeResults.genotypeCommandResult.row;
	var subtypeRow = subtypeRows[0].value;
	var subtypeResult = subtypeRow[1]

	//glue.log("INFO", "subtype RESULT WAS ", subtypeResult);
	if (subtypeResult) {

		//var genoResultElements = subtypeResult.split('_');
		//var subtype = genoResultElements[1];
		//var clade = genoResultElements[1];
		var subtype = subtypeResult.replace("AL_IBV_SEG4_", "");
		glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
		
			glue.command(["set", "field", "rec_subtype", subtype]);
		});
	
	}
	processed++;

	if(processed % 10 == 0) {
		glue.logInfo("Classified "+processed+" segment 4 sequences. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}

});	
