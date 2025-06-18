// Infer Missing gb_country / m49_country from Isolate Name

var placeToCountry = {

  "Sapporo": { country: "Japan", m49: "JPN" },
  "Kanagawa": { country: "Japan", m49: "JPN" },
  "Kyoto": { country: "Japan", m49: "JPN" },
  "Nara": { country: "Japan", m49: "JPN" },
  "Aomori": { country: "Japan", m49: "JPN" },
  "Yamagata": { country: "Japan", m49: "JPN" },
  "Aichi": { country: "Japan", m49: "JPN" },
  "Sendai": { country: "Japan", m49: "JPN" },
  "Mie": { country: "Japan", m49: "JPN" },
  "Osaka": { country: "Japan", m49: "JPN" },
  "Hiroshima": { country: "Japan", m49: "JPN" },
  "Fukuoka": { country: "Japan", m49: "JPN" },
  "Saitama": { country: "Saitama", m49: "JPN" },
  "Miyagi": { country: "Saitama", m49: "JPN" },
  "Shizuoka": { country: "Saitama", m49: "JPN" },
  "Hyogo": { country: "Saitama", m49: "JPN" },

  "Leyte": { country: "Philippines", m49: "PHL" },
  "Biliran": { country: "Philippines", m49: "PHL" },
  "Palawan": { country: "Philippines", m49: "PHL" },

  "Johannesburg": { country: "South Africa", m49: "ZAF" },

  "Sao Paulo": { country: "Brazil", m49: "BRA" },

  "Paris": { country: "France", m49: "FRA" },
  "Berlin": { country: "Germany", m49: "DEU" },
  "Greece": { country: "Greece", m49: "GRC" },
  "England": { country: "United Kingdom", m49: "GBR" },

  "Georgia": { country: "USA", m49: "USA" },
  "Minnesota": { country: "USA", m49: "USA" },
  "Kansas": { country: "USA", m49: "USA" },
  "New Jersey": { country: "USA", m49: "USA" },
  "California": { country: "USA", m49: "USA" },
  "Ann Arbor": { country: "USA", m49: "USA" },
  "Mississippi": { country: "USA", m49: "USA" },

  // Extend this list as needed
};


var result = glue.command([
  "list", "sequence", "sequenceID", "isolate",
  "-w", "source.name = 'icv-ncbi-nuccore' and gb_country = null"
]);
if (result.listResult && result.listResult.row) {
  var seqIDIdx = result.listResult.column.indexOf("sequenceID");
  var isolateIdx = result.listResult.column.indexOf("isolate");
  result.listResult.row.forEach(function(row) {
    var sequenceID = row.value[seqIDIdx];
    var isolate = row.value[isolateIdx];
	if (isolate != null) {
	  var match = isolate.match(/C\/([^\/]+)/);
	  if (match) {
		var place = match[1];
		if (placeToCountry[place]) {
		  glue.inMode("sequence/icv-ncbi-nuccore/" + sequenceID, function() {
		  
			glue.command(["set", "field", "gb_place_sampled", place]);
			glue.command(["set", "field", "gb_country", placeToCountry[place].country]);			
			glue.command(["set", "link-target", "m49_country", "iso_alpha3", placeToCountry[place].m49]);
			
		  });
		  glue.log("INFO", "Set geo info for " + sequenceID + " based on place: " + place);
		} else {
		  glue.log("WARNING", "Place not in lookup table: " + place + " for " + sequenceID);
		}
	  } else {
		glue.log("WARNING", "Could not parse place from isolate: " + isolate + " for " + sequenceID);
	  }
	} else {
	  glue.log("WARNING", "Isolate is null for " + sequenceID);
	}
	
  });
} else {
  glue.log("INFO", "No sequences missing gb_country.");
}

