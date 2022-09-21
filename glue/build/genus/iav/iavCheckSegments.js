var segments = ["1", "2", "3", "4", "5", "6", "7", "8" ];

var sequenceMap = {};

_.each(segments, function(segment) {
	var whereClause = "source.name = 'ncbi-ibv-nuccore' and gb_segment = '"+segment+"'";
	
	var segmentSequences = glue.command(["list", "sequence", "-w", whereClause, "source.name", "sequenceID"], {convertTableToObjects: true});
	
	_.each(segmentSequences, function(segmentSequence) {
		glue.log("INFO", "Sequence Seg");
		sequenceMap[segmentSequence["source.name"]+"/"+segmentSequence["sequenceID"]] = {assignedSegment: "S"+segment, categories:[]}
	});

});

_.each(segments, function(segment) {
	var whereClause = "source.name = 'ncbi-ibv-nuccore' and gb_segment = '"+segment+"'";

	var recognitionResults;
	
	glue.inMode("/module/iavSegmentRecogniser", function() {
		recognitionResults = glue.command(["recognise", "sequence", "-w", whereClause], {convertTableToObjects: true});
	});

	_.each(recognitionResults, function(recognitionResult) {
		var sequenceEntry = sequenceMap[recognitionResult.querySequenceId];
		glue.log("INFO", "Sequence", sequenceEntry);
		sequenceEntry.categories.push({category:recognitionResult.categoryId, direction:recognitionResult.direction});
	});

});

_.each(_.pairs(sequenceMap), function(pair) {
	var seqId = pair[0];
	var seqEntry = pair[1];
	if(seqEntry.categories.length != 1 || !_.isEqual(seqEntry.categories[0], {category:seqEntry["assignedSegment"], direction:"FORWARD"})) {
		glue.log("INFO", "Problematic segment assignment", pair);
	}
});