var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];

var isolateObjs;
glue.inMode("module/tabularUtility", function() {
	isolateObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/core/flu-refseqs.isolates.tsv"]));
});
//glue.logInfo("isolateObjs", isolateObjs);

var isolateIDtoisolateObjs = {};
var i = 2;

_.each(isolateObjs, function(isolateObj) {

	glue.logInfo("row", i);
	//glue.logInfo("isolateObj", isolateObj);

	i++;

	var isolatePK = isolateObjToisolatePK(isolateObj);
	glue.logInfo("isolate public key", isolatePK);

	var virus = isolateObj["virus_name"];
	
	glue.command(["create", "custom-table-row", "isolate", isolatePK]);

	_.each(segments, function(segment) {

		// Skip segent 8 for ICV and IDV
		if (virus == 'icv' && segment == 8 || virus == 'idv' && segment == 8 ) {
	
			// Do nothing
			return;
		}
		else {

			var key = 'segment' + segment + '_accession';
			var segmentSeqID = isolateObj[key];
		
			var sourceName = virus + '-ncbi-refseqs-seg' + segment;
	
			glue.inMode("custom-table-row/isolate/"+isolatePK, function() {


				//glue.logInfo("Linking sequence", segmentSeqID);
				glue.command(["add", "link-target", "sequence", "sequence/"+sourceName+"/"+segmentSeqID]);

				var metadataFields = [
				                      "isolate_id",
									  "iso_source",
									  "iso_country",
									  "iso_region",
									  "iso_year",
									  "iso_month",
									  "iso_day",
									  "iso_host",
									  "lab_host",
									  "cg_subtype",
									  "gb_subtype"];

				_.each(metadataFields, function(metadataField) {

					var value = handleNull(isolateObj[metadataField]);
					if (value) {
				
						//glue.logInfo("isolateObj", isolateObj)
						//glue.logInfo("field", metadataField);
						//glue.logInfo("value", value);
						glue.command(["set", "field", metadataField, value]);

					}

				});
		
				glue.command(["set", "field", 'complete_genome', 'TRUE']);

			});

		}
	
	});

});




// Subroutines

function isolateObjToisolatePK(isolateObj) {
	
	var isolatePK = isolateObj["isolate_id"];
	return isolatePK.replace(/\//g, '|');
	
}

function handleNull(input) {

	if(input == '-') {
		return null;
	}
	return input;
}

