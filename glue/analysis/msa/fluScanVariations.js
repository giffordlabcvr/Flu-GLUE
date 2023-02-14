var viruses = [ 'IAV', 'IBV', 'ICV', 'IDV' ];
var alignmentSets = loadAlignments();
//glue.log("INFO", "RESULT WAS ", alignmentSets);

// Iterate through viruses
_.each(viruses,function(virus) {

	// list variation -w "name like '%phos%' and featureLoc.referenceSequence.name like '%_IAV_%'"
	// list variations
	var where = "name like '%phos%' and featureLoc.referenceSequence.name like '%_" + virus + "_%'";
	var variationResultMap = glue.command(["list", "variation","-w", where]);
	var variationListObj = variationResultMap["listResult"];
	var variationList = variationListObj["row"];	
	//glue.log("INFO", "RESULT WAS ", variationList);

	// Iterate through variations
	_.each(variationList,function(variationObj) {

		var valueArray = variationObj["value"];
		var refSeqName    = valueArray[0];
		var featureName   = valueArray[1];
		var variationName = valueArray[2];
	
		if (alignmentSets[refSeqName]) {
	
			glue.logInfo("Processing variation: "+variationName);

			refseqAlignmentList = alignmentSets[refSeqName];
			
			// Iterate through  alignments 
			for(var i = 0; i < refseqAlignmentList.length; i++) {
			
				var alignmentName = refseqAlignmentList[i]
				
				glue.inMode("alignment/"+alignmentName, function() {
				
					scanResult = glue.tableToObjects(glue.command(["variation", "member", "scan", "-r", refSeqName, "-f", featureName, "-v", variationName]));
					//glue.log("INFO", "load result was:", scanResult);
					
				});

				// Iterate through results
				_.each(scanResult, function(resultObj) {

					//glue.log("INFO", "load result was:", resultObj);

					var sequenceID   = resultObj["sequenceID"];
					var suffCoverage = resultObj["sourceName"];
					var suffCoverage = resultObj["sufficientCoverage"];
					var varIsPresent = resultObj["present"];
					//glue.logInfo("TEST: "+varIsPresent+" in "+sequenceID);
					
					if (varIsPresent) {
					
						// Update the variation field
						glue.logInfo("Variation: "+variationName+" is present in "+sequenceID);
	
					}
					else {
						
						glue.logInfo("Variation: "+variationName+" is NOT present in "+sequenceID);
					
					
					}
	
				});
			
			}
			
		}
		else {
			
			glue.logInfo("Skipping variation"+variationName+" - it's reference "+refSeqName+" does not constrain any alignments");
		
		}
	
	});

});


// Subroutines
function loadAlignments() {

	var alignmentSets = {};
	
	
	// list alignments
	var where = "name not like '%_ROOT_%'";
	var alignmentResultMap = glue.command(["list", "alignment","-w", where]);
	var alignmentListObj = alignmentResultMap["listResult"];
	var alignmentList = alignmentListObj["row"];
	//glue.log("INFO", "RESULT WAS ", alignmentList);

	_.each(alignmentList,function(alignmentObj){

		var valueArray    = alignmentObj["value"];
		var alignmentName = valueArray[0];
		var refSeqName    = valueArray[2];	
		glue.log("INFO", "RESULT WAS ", refSeqName);

		if (alignmentSets[refSeqName]) {
		
			refseqAlignmentList = alignmentSets[refSeqName];
			glue.log("INFO", "RESULT WAS ", refseqAlignmentList);
			refseqAlignmentList.push(alignmentName);
			
		}
		else {
			var refseqAlignmentList = [];
			refseqAlignmentList.push(alignmentName);					
			alignmentSets[refSeqName] = refseqAlignmentList;		
		}
		
	});

	return alignmentSets;

}
