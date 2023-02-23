var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];

var strainObjs;

glue.inMode("module/tabularUtility", function() {
	strainObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/genus/iav/iav-strains.complete.out.tsv"]));
});
//glue.logInfo("strainObjs", strainObjs);

var strainIDtoStrainObjs = {};
var i = 2;

_.each(strainObjs, function(strainObj) {

	glue.logInfo("row", i);
	//glue.logInfo("StrainObj", strainObj); die;
	
	i++;
	
	var strainPK = strainObjToStrainPK(strainObj);

 	glue.command(["create", "custom-table-row", "isolate", strainPK]);


	_.each(segments, function(segment) {

		var key = 'segment' + segment + '_accession';
		var segmentSeqID = strainObj[key];		
		var sourceName = 'iav-ncbi-curated-segment-' + segment;
		
	
		glue.inMode("custom-table-row/isolate/"+strainPK, function() {


			glue.logInfo("Linking sequence", segmentSeqID);
			glue.command(["add", "link-target", "sequence", "sequence/"+sourceName+"/"+segmentSeqID]);

			var metadataFields = ["strain",
								  "isolate", 
								  "iso_source",
								  "iso_country",
								  "iso_place_name",
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
					
					glue.logInfo("strainObj", strainObj)
					glue.logInfo("field", metadataField);
					glue.logInfo("value", value);
					glue.command(["set", "field", metadataField, value]);

				}

			});

			glue.command(["set", "field", 'complete_genome', 'TRUE']);

		});
		
			
		// Set the recogniser segment field
		glue.inMode("sequence/"+sourceName+"/"+segmentSeqID, function() {

			glue.command(["set", "field", "rec_segment", segment]);
			
		});

		var completeGenomeSubtype = strainObj["cg_subtype"];	
		var genomeSubtypeArray = completeGenomeSubtype.split("|");
		//glue.logInfo("genomeSubtypeArray", genomeSubtypeArray); die;

		glue.inMode("sequence/"+sourceName+"/"+segmentSeqID, function() {

			var index = segment - 1;
			var recSubtype = genomeSubtypeArray[index];
			glue.command(["set", "field", "rec_subtype", recSubtype]);				

		});
		
		
	});

});


// Subroutines

function strainObjToStrainPK(strainObj) {
	
	var strainPK = strainObj["strain_id"];
	return strainPK.replace(/\//g, '|');
	
}

function handleNull(input) {

	if(input == '-') {
		return null;
	}
	return input;
}

