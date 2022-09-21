// Script to calculate nucleotide composition in reference sequences
aaCompositionResults = {}

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

// iterate through reference list
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
// Create column headers
// Reference Feature A C D E F etc
var headerRow1 = [  'ReferenceName', 'FeatureName' ];
var headerRowAas = ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P' ];
var outputArray = [];

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


        var rowStart = [  referenceName, featureName, length ];
        var rowFirstPart = rowStart.join("\t");
        var aaResults = [];
        _.each(headerRowAas, function(aa) {
        
        	var aaCount = aaObj[aa];
        	var aaFreq = aaCount / length;
        	var aaFormatedFreq = aaFreq.toFixed(2);
        	//glue.log("INFO", "  Amino acid '"+aa+" frequency = ("+aaCount+" / "+length+") "+aaFormatedFreq);
        	aaResults.push(aaCount);
         	aaResults.push(aaFormatedFreq);

        });

        var rowSecondPart = aaResults.join("\t");
        var completeRow = [ rowFirstPart, rowSecondPart, "\n" ]
        outputArray.push(completeRow);
        glue.log("INFO", "ROW "+completeRow);
	
	 });

 });


    
// Export results








































