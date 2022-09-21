
var allResults = [];

glue.inMode("module/iavSegmentRecogniser", function(){

	 var recogniserResult = glue.command(["recognise","sequence","-a"]);
     recogniserResult = recogniserResult["blastSequenceRecogniserResult"];
	 //glue.log("INFO", "Segment recogniser result was:", recogniserResult);
    

     // Iterate through rows updating recogniser_segment field
     
     
     var tableRows = recogniserResult["row"];
	 //glue.log("INFO", "Segment recogniser result was:", tableRows);


	 _.each(tableRows, function(rowObj)  {

		var valueObj = rowObj["value"];
		var querySequenceId = valueObj[0];
		var recogniserSeg = valueObj[1];
		var direction = valueObj[2];

		var idElements = querySequenceId.split('/');
		var sourceName = idElements[0];
		var sequenceID = idElements[1];
	    //glue.log("INFO", "Got segment '"+recogniserSeg+"' for sequence: '"+sequenceID+"'");
		
		var result = {};
		result["sourceName"] = sourceName;
		result["sequenceID"] = sequenceID;
		result["recogniserSeg"] = recogniserSeg;

        allResults.push(result);

 	 
	 });

});

//glue.log("INFO", "Segment recogniser result was:", allResults);
    

_.each(allResults, function(resultObj)  {

	 // update the sequence table with the results
	 
	 var sourceName = resultObj["sourceName"];	 
	 var sequenceID = resultObj["sequenceID"];	 
	 var recogniserSeg = resultObj["recogniserSeg"];
	 	 
	 glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
	 
		 
	 	 if (recogniserSeg) {
		 	glue.command(["set", "field", "recogniser_segment", recogniserSeg]);
 	 	 }

 	 
	 });
	 
});
