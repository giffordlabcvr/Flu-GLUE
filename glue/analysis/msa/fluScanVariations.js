// project flu list variation -w "name like '%phos%'"

// list features
var codingFeatures = {};
var resultMap = glue.command(["list", "feature","-w", "featureMetatags.name = 'CODES_AMINO_ACIDS'"]);
var featureList = resultMap["listResult"];
var codingFeatureList = featureList["row"];
_.each(codingFeatureList,function(featureObj){

	//glue.log("INFO", "RESULT WAS ", featureObj);

	var valueArray = featureObj["value"];
	var codingFeatureName = valueArray[0];
	//glue.log("INFO", "NAME WAS ", featureName)
	codingFeatures[codingFeatureName] = featureObj;

	// project flu list variation -w "name like '%phos%' and featureLoc.feature.name = 'NP' and featureLoc.referenceSequence.name like '%IAV_MASTER%'"
	var variations;

	// Scan alignments for variations
	var count = 0;
	for(var i = 0; i < variations.length; i++) {


		count++;
	
		var variationName = variations[i]
		glue.logInfo("Processing variation: "+variationName);

		// var alignmentName;


		var scanResult;
		glue.inMode("alignment/"+alignmentName, function() {
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
		
					glue.log("INFO", "Variation name result was:"+variationName+" sequence: "+sequenceID+" first ref codon: "+firstRefCodon);
					
					//glue.command(["set", "field", variationName, firstRefCodon]);	
			
			
				});
			
		
			}

		
		});


	}		


});
