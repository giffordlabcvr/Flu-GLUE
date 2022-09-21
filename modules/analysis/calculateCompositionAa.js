function calculateCompositionAa() {

	// Script to calculate amino acid composition in reference sequences

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

	
	});	
	//glue.log("INFO", "RESULT WAS ", codingFeatures);

	// get list of reference sequences from GLUE
	var referencesResult = glue.command(["list","reference"]);
	//glue.log("INFO", "RESULT WAS ", referencesResult);

	var listResult = referencesResult["listResult"];
	var referencesList = listResult["row"];
	//glue.log("INFO", "RESULT WAS ", referencesList);


	// iterate through reference list and get AA composition of each coding feature
	aaCompositionResults = {}
	_.each(referencesList, function(refObj) {

		//glue.log("INFO", "RESULT WAS ", refObj);
	
		var refseqResults = {};
		var referenceProperties = refObj["value"];
		var referenceName = referenceProperties[0];
	
		//glue.log("INFO", "Reference name result was:", referenceName);

		// list all features annotated in this reference 
		// GLUE COMMAND: reference [referenceName] list feature-location
		glue.inMode("/reference/"+referenceName, function() {

			var featuresResult = glue.tableToObjects(glue.command(["list", "feature-location"]));
			//glue.log("INFO", "RESULT WAS ", featuresResult);
			 
			// iterate through features
			_.each(featuresResult, function(featureObj) {

			   //glue.log("INFO", "RESULT WAS ", featureObj);
			   var featureResults = {};
		   
			   // get amino acid sequence
			   var featureName = featureObj["feature.name"];
			   //glue.log("INFO", "Feature name result was:", featureName);
		   
			   // get amino acid sequence of the features
			   if (codingFeatures[featureName]) {

				   glue.inMode("/feature-location/"+featureName, function() {
						  
					   var aaSequenceResult = glue.tableToObjects(glue.command(["amino-acid"]));

					   //glue.log("INFO", "Amino acid seq result was:", aaSequenceResult);
	 
					   // iterate through and get amino acid composition of feature
					   var length = 0;
					   _.each(aaSequenceResult, function(codonObj) {
 
						   var aa = codonObj["aminoAcid"];
					   
						   length += 1;

						   //glue.log("INFO", "Amino acid residue result was:", aa);
			
						   if (featureResults[aa]) {

							   featureResults[aa] += 1;
			
						   }
						   else {

							   featureResults[aa] = 1;
						   }

					   });
				   
					   featureResults["length"] = length;
		   
				   });
			   
				  // store feature result
				  refseqResults[featureName] = featureResults;
			
				}
				   
			});   

			// store reference result
			aaCompositionResults[referenceName] = refseqResults;
		
		});   
	
	});

	//glue.log("INFO", "FINAL RESULT WAS ", aaCompositionResults);

	// Transform the data 
	var outputArray = [];
	var headerRowAas = ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P' ];

	// Iterate through reference sequences
	_.each(_.keys(aaCompositionResults), function(referenceName) {

		//glue.log("INFO", "Got reference name '"+referenceName);
	
		// Iterate through reference sequence features
		var featuresObj = aaCompositionResults[referenceName];
   
		_.each(_.keys(featuresObj), function(featureName) {

		   // Write values for each amino (count + ratio) 
		   var aaObj = featuresObj[featureName];
		   var length = aaObj["length"];
		   //glue.log("INFO", "Got reference name '"+referenceName+" and feature "+featureName+" length = "+length);
		   
		   
		   var aaCountResults = {};
		   var aaFreqResults = {};
		   _.each(headerRowAas, function(aa) {

			   var aaFreq;
			   var aaCount = aaObj[aa];
			   if (aaCount) {
					var aaFreq = (aaCount / length) * 100;
					var aaFormatedFreq = aaFreq.toFixed(2);
					//glue.log("INFO", "  Amino acid '"+aa+" frequency = ("+aaCount+" / "+length+") "+aaFormatedFreq);
			   }
			   else {
				    aaCount = '0';
				    aaFormatedFreq = '0';		  
			   }
			  
			  aaCountResults[aa] = aaCount;
			  aaFreqResults[aa] = aaFormatedFreq;

		   });

		   // add results to array to be returned to GLUE
		   outputArray.push({
		
			   referenceName: referenceName,
			   featureName: featureName,
			   seqLength: length,
			   A: aaFreqResults["A"],
			   C: aaFreqResults["C"],
			   D: aaFreqResults["D"],
			   E: aaFreqResults["E"],
			   F: aaFreqResults["F"],
			   G: aaFreqResults["G"],
			   H: aaFreqResults["H"],
			   I: aaFreqResults["I"],
			   K: aaFreqResults["K"],
			   L: aaFreqResults["L"],
			   M: aaFreqResults["M"],
			   N: aaFreqResults["N"],
			   P: aaFreqResults["P"],
			   Q: aaFreqResults["Q"],
			   R: aaFreqResults["R"],
			   S: aaFreqResults["S"],
			   T: aaFreqResults["T"],
			   V: aaFreqResults["V"],
			   W: aaFreqResults["W"],
			   Y: aaFreqResults["Y"]

			   //A: aaFreqResults["A"]+" ("+aaCountResults["A"]+")",
			   //C: aaFreqResults["C"]+" ("+aaCountResults["C"]+")",
			   //D: aaFreqResults["D"]+" ("+aaCountResults["D"]+")",
			   //E: aaFreqResults["E"]+" ("+aaCountResults["E"]+")",
			   //F: aaFreqResults["F"]+" ("+aaCountResults["F"]+")",
			   //G: aaFreqResults["G"]+" ("+aaCountResults["G"]+")",
			   //H: aaFreqResults["H"]+" ("+aaCountResults["H"]+")",
			   //I: aaFreqResults["I"]+" ("+aaCountResults["I"]+")",
			   //K: aaFreqResults["K"]+" ("+aaCountResults["K"]+")",
			   //L: aaFreqResults["L"]+" ("+aaCountResults["L"]+")",
			   //M: aaFreqResults["M"]+" ("+aaCountResults["M"]+")",
			   //N: aaFreqResults["N"]+" ("+aaCountResults["N"]+")",
			   //P: aaFreqResults["P"]+" ("+aaCountResults["P"]+")",
			   //Q: aaFreqResults["Q"]+" ("+aaCountResults["Q"]+")",
			   //R: aaFreqResults["R"]+" ("+aaCountResults["R"]+")",
			   //S: aaFreqResults["S"]+" ("+aaCountResults["S"]+")",
			   //T: aaFreqResults["T"]+" ("+aaCountResults["T"]+")",
			   //V: aaFreqResults["V"]+" ("+aaCountResults["V"]+")",
			   //W: aaFreqResults["W"]+" ("+aaCountResults["W"]+")",
			   //Y: aaFreqResults["Y"]+" ("+aaCountResults["Y"]+")"



		   });

   
		});


	});
	 
	return outputArray;
}


