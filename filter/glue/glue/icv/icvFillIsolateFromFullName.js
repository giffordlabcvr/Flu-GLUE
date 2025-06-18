// === Fill missing 'isolate' from 'full_name' field ===

var result = glue.command([
  "list", "sequence", "sequenceID", "full_name", "isolate",
  "-w", "source.name = 'icv-ncbi-nuccore' and isolate = null"
]);
if (result.listResult && result.listResult.row) {
  var seqIDIdx = result.listResult.column.indexOf("sequenceID");
  var nameIdx = result.listResult.column.indexOf("full_name");
  result.listResult.row.forEach(function(row) {
    var sequenceID = row.value[seqIDIdx];
    var fullName = row.value[nameIdx];
    var match = fullName.match(/\((C\/[^)]+)\)/);  // Adjust for species as needed
    if (match) {
      var inferredIsolate = match[1];
      glue.inMode("sequence/icv-ncbi-nuccore/" + sequenceID, function() {
        glue.command(["set", "field", "isolate", inferredIsolate]);
      });
      glue.log("INFO", "Set isolate for " + sequenceID + ": " + inferredIsolate);
    } else {
      glue.log("WARNING", "Could not parse isolate from full_name for " + sequenceID);
    }
  });
} else {
  glue.log("INFO", "No sequences missing isolate.");
}

