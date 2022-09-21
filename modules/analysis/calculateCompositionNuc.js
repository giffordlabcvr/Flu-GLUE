// Calculate nucleotide composition in reference sequences
function calculateCompositionNuc() {

	// Script to calculate nucleotide composition in reference sequences
	var nucCompositionResults = {};
	var lengthResults = {};

	// export reference sequences from GLUE
	var outputArray = [];
	glue.inMode("module/fastaExporter", function(){

		var sequences = glue.command(["export","-w","source.name like '%refseq%'","-p"]);
		var list = sequences.nucleotideFasta.sequences;
	
		_.each(list, function(seq)  {

			//loop through each sequence in the alignment
			seqNucComposition = {}
		
			var sequence   = seq.sequence;
			var sequenceId = seq.id;
			//glue.log("INFO", "ID result was:", sequenceId);

			//loop through each position in the current sequence
			for (var i=0; i < seq.sequence.length; i++) {
				var base = seq.sequence.charAt(i);
			
				if (seqNucComposition[base]) {
					seqNucComposition[base] += 1;		
				}
				else {
					seqNucComposition[base] = 1;
				}	
			}
		
			//glue.log("INFO", "Nucleotide composition result was:", seqNucComposition);
		
			nucCompositionResults[sequenceId] = seqNucComposition;
			lengthResults[sequenceId] = seq.sequence.length;
			
		});

		_.each(_.keys(nucCompositionResults), function(sequenceID) {
	
			var seqResults = nucCompositionResults[sequenceID];
			//glue.log("INFO", "RESULTS FOR '"+sequenceID+"'");

		    var nucCountResults = {};
		    var nucFreqResults = {};
            var length = lengthResults[sequenceID];


			_.each(_.keys(seqResults), function(nucleotide) {

				var count = seqResults[nucleotide];
				var length = lengthResults[sequenceID];		
				var ratio = (count / length) * 100;
				var formatRatio = ratio.toFixed(2);
				//glue.log("INFO", "Got count '"+count+"' for nucleotide character'"+nucleotide+"'");
				//glue.log("INFO", "Got ratio '"+ratio+"' for nucleotide character'"+nucleotide+"'");
				
				nucCountResults[nucleotide] = count;
				nucFreqResults[nucleotide] = formatRatio;
				
		
			});

			// add results to array to be returned to GLUE
			outputArray.push({
		
				sequenceID: sequenceID,
				seqLength: length,
				"A#": nucCountResults["A"],
				"T#": nucCountResults["T"],
				"G#": nucCountResults["G"],
				"C#": nucCountResults["C"],
				"A%": nucFreqResults["A"],
				"T%": nucFreqResults["T"],
				"G%": nucFreqResults["G"],
				"C%": nucFreqResults["C"],
			});

		
		});

	});
 	 
	return outputArray;
}


