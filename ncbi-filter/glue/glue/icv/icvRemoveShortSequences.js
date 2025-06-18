// Remove short sequences

var result = glue.command([
  "list", "sequence", "sequenceID",
  "-w", "length < 300 and source.name = 'icv-ncbi-nuccore'"
]);

if (result.listResult && result.listResult.row) {
  var colIdx = result.listResult.column.indexOf("sequenceID");
  result.listResult.row.forEach(function(row) {
    var sequenceID = row.value[colIdx];
    glue.command(["delete", "sequence", "-w", "sequenceID = '"+sequenceID+"'"]);
    glue.log("INFO", "Deleted short sequence: " + sequenceID);
  });
} else {
  glue.log("INFO", "No short sequences found to delete.");
}
