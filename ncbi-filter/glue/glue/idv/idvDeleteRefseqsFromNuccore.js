// === Delete sequences from a source if they appear in the reference set ===

var sourceName = "idv-ncbi-nuccore";

// Step 1: Collect sequenceIDs from refseq sources
var refResult = glue.command([
  "list", "sequence", "sequenceID",
  "-w", "source.name like 'idv-ncbi-refseqs-%'"
]);

var refIDs = {};
if (refResult.listResult && refResult.listResult.row) {
  var idIdx = refResult.listResult.column.indexOf("sequenceID");
  refResult.listResult.row.forEach(function(row) {
    var seqID = row.value[idIdx];
    refIDs[seqID] = true;
  });
} else {
  glue.log("WARNING", "No sequences found in refseq sources.");
}

// Step 2: Delete duplicates from nuccore
var seqResult = glue.command([
  "list", "sequence", "sequenceID",
  "-w", "source.name = '" + sourceName + "'"
]);

if (seqResult.listResult && seqResult.listResult.row) {
  var idIdx = seqResult.listResult.column.indexOf("sequenceID");
  seqResult.listResult.row.forEach(function(row) {
    var seqID = row.value[idIdx];
    if (refIDs[seqID]) {
      glue.command([
        "delete", "sequence",
        "-w", "source.name = '" + sourceName + "' and sequenceID = '" + seqID + "'"
      ]);
      glue.log("INFO", "Deleted duplicate sequence from nuccore: " + seqID);
    }
  });
} else {
  glue.log("INFO", "No sequences found in source: " + sourceName);
}
