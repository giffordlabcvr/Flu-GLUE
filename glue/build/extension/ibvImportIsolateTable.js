var segments = [1, 2, 3, 4, 5, 6, 7, 8 ];

var isolateObjs;
glue.inMode("module/tabularUtility", function() {
    isolateObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/extension/ibv_nuccore_isolates.tsv"]));
});

var i = 2;

_.each(isolateObjs, function(isolateObj) {

    glue.logInfo("row", i);
    i++;

    var isolatePK = isolateObjToisolatePK(isolateObj);
    glue.logInfo("isolate public key", isolatePK);

    var species = isolateObj["species"];
    var isolatePath = "custom-table-row/isolate/" + isolatePK;

    // Check if isolate row already exists
    var exists = false;
    try {
        glue.inMode(isolatePath, function() {
            exists = true;
        });
    } catch (e) {
        // not found, safe to create
    }

    if (exists) {
        glue.logInfo("Skipping existing isolate", isolatePK);
        return;
    }

    glue.command(["create", "custom-table-row", "isolate", isolatePK]);

    glue.inMode(isolatePath, function() {

        // Set metadata fields (unconditionally)
        var metadataFields = [
            "isolate_id", "species", "origin_type", "sample_type", "iso_source", "iso_country",
            "iso_region", "iso_year", "iso_month", "iso_day", "host", "lab_host",
            "gb_serotype", "rec_serotype", "mlca_serotype", "genome_lineage",
            "segment1_accession", "segment2_accession", "segment3_accession", "segment4_accession",
            "segment5_accession", "segment6_accession", "segment7_accession", "segment8_accession"
        ];

        _.each(metadataFields, function(metadataField) {
            var value = handleNull(isolateObj[metadataField]);

			if (value !== null && value !== undefined) {
				//glue.logInfo("Setting field", metadataField + " = " + value);
				
				// Wrap in quotes if the value contains the pipe character 
				if (/[|\s]/.test(value)) {
					value = '"' + value + '"';
				}
				glue.command(["set", "field", metadataField, value]);				
			}


        });

        // Count found segments
        var foundSegments = 0;
        _.each(segments, function(segment) {
            var key = 'segment' + segment + '_accession';
            var value = isolateObj[key];
            if (value && value !== '-') {
                foundSegments++;
            }
        });

        glue.command(["set", "field", "segment_count", foundSegments.toString()]);
        glue.command(["set", "field", "is_complete", (foundSegments === 8 ? "TRUE" : "FALSE")]);

		// Flag this isolate as 'rank-and-file' (i.e. not represented by the core reference sequences)
		glue.command(["set", "field", "is_reference", "FALSE"]);

    });

});


// Subroutines
function isolateObjToisolatePK(isolateObj) {
    return isolateObj["isolate_id"].replace(/\//g, '|');
}

function handleNull(input) {
    return (input == '-' ? null : input);
}
