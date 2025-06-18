// Remove short sequences

var result = glue.command([
  "list", "sequence", "sequenceID",
  "-w", "patent_related = 'true'"
]);

if (result.listResult && result.listResult.row) {
  var colIdx = result.listResult.column.indexOf("sequenceID");
  result.listResult.row.forEach(function(row) {
    var sequenceID = row.value[colIdx];
    glue.command(["delete", "sequence", "-w", "sequenceID = '"+sequenceID+"'"]);
    glue.log("INFO", "Deleted patent sequence: " + sequenceID);
  });
} else {
  glue.log("INFO", "No patent sequences found to delete.");
}
