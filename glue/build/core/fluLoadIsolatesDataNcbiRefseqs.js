var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];

var strainObjs;
glue.inMode("module/tabularUtility", function() {
	strainObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/core/flu-refseqs.isolates.tsv"]));
});
//glue.logInfo("strainObjs", strainObjs);

var strainIDtoStrainObjs = {};
var i = 2;

_.each(strainObjs, function(strainObj) {

	glue.logInfo("row", i);
	//glue.logInfo("StrainObj", strainObj);

	i++;

	var strainPK = strainObjToStrainPK(strainObj);
	var completeGenomeSubtype = strainObj["cg_subtype"];
	var virus = strainObj["virus_name"];
	
	var genomeSubtypeArray = completeGenomeSubtype.split("|");
	//glue.logInfo("genomeSubtypeArray", genomeSubtypeArray); die;

	glue.command(["create", "custom-table-row", "isolate", strainPK]);


	_.each(segments, function(segment) {

		// Skip segent 8 for ICV and IDV
		if (virus == 'icv' && segment == 8) {
	
		
		}
		else if (virus == 'idv' && segment == 8) {
	
			// Do nothing
		}
		else {

			var key = 'segment' + segment + '_accession';
			var segmentSeqID = strainObj[key];
		
			var sourceName = virus + '-ncbi-refseqs-seg' + segment;
	
			glue.inMode("custom-table-row/isolate/"+strainPK, function() {


				//glue.logInfo("Linking sequence", segmentSeqID);
				glue.command(["add", "link-target", "sequence", "sequence/"+sourceName+"/"+segmentSeqID]);

				var metadataFields = [
									  "iso_source",
									  "iso_country",
									  "iso_region",
									  "iso_year",
									  "iso_month",
									  "iso_day",
									  "iso_host",
									  "lab_host",
									  "cg_subtype",
									  "hn_subtype"];

				_.each(metadataFields, function(metadataField) {

					var value = handleNull(strainObj[metadataField]);
					if (value) {
				
						//glue.logInfo("strainObj", strainObj)
						//glue.logInfo("field", metadataField);
						//glue.logInfo("value", value);
						glue.command(["set", "field", metadataField, value]);

					}

				});
		
				glue.command(["set", "field", 'complete_genome', 'TRUE']);

			});


			// Set the recogniser segment field
			glue.inMode("sequence/"+sourceName+"/"+segmentSeqID, function() {
			
				glue.command(["set", "field", "rec_segment", segment]);
						
			});

			// Set the recogniser subtype field for IAV
			if (virus == 'iav') {
				glue.inMode("sequence/"+sourceName+"/"+segmentSeqID, function() {
			
					var index = segment - 1;
					var recSubtype = genomeSubtypeArray[index];
					glue.command(["set", "field", "rec_subtype", recSubtype]);				

				});
				
			}

		}
	
	});

});




// Subroutines

function strainObjToStrainPK(strainObj) {
	
	var strainPK = strainObj["strain_id"];
	glue.logInfo("strain public key", strainPK);

	return strainPK.replace(/\//g, '_');
	
}

function handleNull(input) {

	if(input == '-') {
		return null;
	}
	return input;
}

