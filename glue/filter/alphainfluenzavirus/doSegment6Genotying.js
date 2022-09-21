
var allResults = [];

glue.inMode("module/iavSubtypeRecogniserSegment6", function(){

	 var recogniserResult = glue.command(["recognise","sequence","-w", "recogniser_segment = 6"]);
     recogniserResult = recogniserResult["blastSequenceRecogniserResult"];
	 //glue.log("INFO", "Segment recogniser result was:", recogniserResult);
    

     // Iterate through rows updating recogniser_segment field
     
     
     var tableRows = recogniserResult["row"];
	 //glue.log("INFO", "Segment recogniser result was:", tableRows);


	 _.each(tableRows, function(rowObj)  {

		var valueObj = rowObj["value"];
		var querySequenceId = valueObj[0];
		var recogniserLineage = valueObj[1];
		var direction = valueObj[2];

		var idElements = querySequenceId.split('/');
		var sourceName = idElements[0];
		var sequenceID = idElements[1];
	    glue.log("INFO", "Got subtype '"+recogniserLineage+"' for sequence: '"+sequenceID+"'");
		
		var result = {};
		result["sourceName"] = sourceName;
		result["sequenceID"] = sequenceID;
		result["recogniserLineage"] = recogniserLineage;

        allResults.push(result);

 	 
	 });

});

//glue.log("INFO", "Segment recogniser result was:", allResults);
_.each(allResults, function(resultObj)  {

	 // update the sequence table with the results
	 
	 var sourceName = resultObj["sourceName"];	 
	 var sequenceID = resultObj["sequenceID"];	 
	 var recogniserLineage = resultObj["recogniserLineage"];
	 	 
	 glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
	 
	 	 if (recogniserLineage) {
		 	glue.command(["set", "field", "recogniser_lineage", recogniserLineage]);
 	 	 }
 	 
	 });
	 
});
