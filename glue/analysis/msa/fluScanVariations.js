
// Scan alignments for variations
var variations = [ "n-linked-glycosylation-o1", "n-linked-glycosylation-o2", "n-linked-glycosylation-o3", "n-linked-glycosylation-o4", "n-linked-glycosylation-o5", "n-linked-glycosylation-o6" ];

var count = 0;
for(var i = 0; i < variations.length; i++) {


	count++;

	var variationName = variations[i]
	glue.logInfo("Processing variation: "+variationName);

	var scanResult;
	glue.inMode("alignment/AL_IFNL_MAMMAL_B", function() {
		scanResult = glue.tableToObjects(glue.command(["variation", "member", "scan", "-r", "REF_IFNL_Mammal_b_MASTER", "-f", "orf", "-v", variationName, "-t"]));
		//glue.log("INFO", "load result was:", scanResult);
	});

	// Iterate through results
	_.each(scanResult, function(resultObj) {
		
		var sequenceID = resultObj["sequenceID"];
		var sourceName   = resultObj["sourceName"];
		var firstRefCodon   = resultObj["firstRefCodon"];

		if (firstRefCodon) {
	
			// Update the variation field
			glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
	
 				var variationName = 'glycosylation_o' + count;
		
				glue.log("INFO", "Vriation name result was:"+variationName+" sequence: "+sequenceID+" first ref codon: "+firstRefCodon);
				glue.command(["set", "field", variationName, firstRefCodon]);	
			
			
			});
			
		
		}

		
	});


}		



