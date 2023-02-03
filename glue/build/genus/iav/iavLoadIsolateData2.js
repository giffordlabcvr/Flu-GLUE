// First load the strain data
var strainDatabyStrainId = {};
load_strain_data(strainDatabyStrainId);
glue.log("INFO", "By strain hash", strainDatabyStrainId);


var isolate_fields = [ 'iso_host', 'iso_country', 'iso_place_name', 
                       'rec_cg_subtype', 'rec_hn_subtype', 'iso_day' ];
var segments = [ 1, 2, 3, 4, 5, 6, 7, 8 ];


_.each(segments, function(segment) {

    var sourceName = 'iav-ncbi-curated-segment-' + segment;
    
	// list the sequences in source
	var listSeqResult = glue.command(["list", "sequence", "-w", "source.name = '"+sourceName+"'"]);

	// extract from the result a list of sequence IDs.
	var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");
	glue.log("INFO", "Source name", sourceName);
	
	// for each sequence ID
	_.each(seqIds, function(seqId) {

		// create an object in the custom table which uses the sequence ID as the row ID.
		glue.command(["create", "custom-table-row", "isolate_data", seqId]);
	
		// associate the corresponding sequence with this object.
		glue.inMode("sequence/"+sourceName+"/"+seqId, function() {
		
			glue.command(["set", "link-target", "isolate_data", "custom-table-row/isolate_data/"+seqId]);
			glue.log("INFO", "Seqeunce", seqId);

			glue.command(["set", "field", "name", 'IAV']);
			
		});

		_.each(isolate_fields, function(isolate_field) {
			
			glue.log("INFO", "FIELD", isolate_field);
			var strain_isolation_data = strainDatabyStrainId[seqId];		
			glue.log("INFO", "strain_isolation_data", strain_isolation_data);
			var iso_field_value = strain_isolation_data[isolate_field];
			if (iso_field_value) {
				glue.command(["custom-table-row", "isolate_data", seqId, "set", "field", isolate_field, iso_field_value]);
			    die;			
			}
			else {
				glue.log("INFO", "NO VALUE FOT TABLE FIELD", isolate_field);
				die;
			
			}
			
		});

	});

});



_.each(segments, function(segment) {

    var sourceName = 'iav-ncbi-refseqs-seg' + segment;
    
	// list the sequences in source
	var listSeqResult = glue.command(["list", "sequence", "-w", "source.name = '"+sourceName+"'"]);

	// extract from the result a list of sequence IDs.
	var seqIds = glue.getTableColumn(listSeqResult, "sequenceID");
	//glue.log("INFO", "Source name", sourceName);

	// for each sequence ID
	_.each(seqIds, function(seqId) {

		// create an object in the custom table which uses the sequence ID as the row ID.
		glue.command(["create", "custom-table-row", "isolate_data", seqId]);
	
		// associate the corresponding sequence with this object.
		glue.inMode("sequence/"+sourceName+"/"+seqId, function() {
		
			glue.command(["set", "link-target", "isolate_data", "custom-table-row/isolate_data/"+seqId]);
			glue.command(["set", "field", "name", 'IAV']);	
						
		});

		// update isolate data in the DB
		glue.inMode("custom-table-row/isolate_data/"+seqId, function() {
		

			_.each(isolate_fields, function(isolate_field) {
				
				var iso_field_value = strainDatabyStrainId[isolate_field];
				if (iso_field_value) {
				
					glue.command(["set", "field", isolate_field, iso_field_value]);
				
				}
    			
    		});

		});

	});

});


// Subroutine - load strain data to a hash of hashes (field=value) with sequenceIDs as keys
function load_strain_data(strainDatabyStrainId) {


	var strainData = glue.command(["file-util", "load-string", "tabular/genus/iav/strains.complete.out.tsv"]).fileUtilLoadStringResult.loadedString;

	//glue.log("INFO", "Strain data string", strainData);
	var myArray = strainData.split("\n");
	//glue.log("INFO", "Strain array", myArray);


	var i = 0;
	var column_headers = {};
	_.each(myArray, function(array_string) {

		//glue.log("INFO", "Strain data string", array_string);
		i++;
	
		var myLineDataArray = [];
		myLineDataArray = array_string.split("\t");
		//glue.log("INFO", "Strain data string", myLineDataArray);
		//glue.log("INFO", "poove1", i); 
		if (i == '1') {
		
			var j = 0;
			_.each(myLineDataArray, function(field_name) {
	 
				j++;
				//glue.log("INFO", "Header column", field_name);
				column_headers[j] = field_name;
			});
		
		}  
		else {

			var k = 0;
			var seqID;
			_.each(myLineDataArray, function(value) {


			
				k++;
				var field_name = column_headers[k];

				//glue.log("INFO", "\tfield", field_name);
				//glue.log("INFO", "\tvalue", value);

				if (k == '1') {
			
					seqID = value;
					var seqData = {};
					strainDatabyStrainId[seqID] = seqData;
				}
				else {
					seqData = strainDatabyStrainId[seqID];
					seqData[field_name] = value;			
				}
			
			});
	
		}

	});
	
	
}
