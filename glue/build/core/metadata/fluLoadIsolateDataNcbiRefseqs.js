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

	var species = isolateObj["species"];
	
	glue.command(["create", "custom-table-row", "isolate", isolatePK]);

	var foundSegments = 0;
	var expectedSegments = (species == 'ICV' || species == 'IDV') ? 7 : 8;

	glue.inMode("custom-table-row/isolate/" + isolatePK, function() {

		_.each(segments, function(segment) {
		
			if ((species == 'ICV' || species == 'IDV') && segment == 8) return;

			var key = 'segment' + segment + '_accession';
			var segmentSeqID = isolateObj[key];
			var sourceName = species + '-ncbi-refseqs-seg' + segment;

			if (segmentSeqID) {
				foundSegments++;
				glue.command(["add", "link-target", "sequence", "sequence/" + sourceName + "/" + segmentSeqID]);

				var metadataFields = [
					"isolate_id", "species", "origin_type", "sample_type", "iso_source", "iso_country",
					"iso_region", "iso_year", "iso_month", "iso_day", "host", "lab_host",
					"gb_serotype", "rec_serotype", "mlca_serotype", "genome_lineage"
				];

				_.each(metadataFields, function(metadataField) {
					var value = handleNull(isolateObj[metadataField]);
					if (value) {
						glue.command(["set", "field", metadataField, value]);
					}
				});
			}
		});

		// Set segment count and completeness
		glue.command(["set", "field", "segment_count", foundSegments.toString()]);
		if (foundSegments == expectedSegments) {
			glue.command(["set", "field", "is_complete", "TRUE"]);
		} else {
			glue.command(["set", "field", "is_complete", "FALSE"]);
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

