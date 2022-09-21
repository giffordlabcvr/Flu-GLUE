// Calculate dinucleotide composition in reference sequences
function calculateCompositionDinuc() {

	var dinucCompositionResults = {};
	var lengthResults = {};
    var outputArray = [];
	
	// export reference sequences from GLUE
	glue.inMode("module/fastaExporter", function(){

		var sequences = glue.command(["export","-w","source.name like '%refseqs%'","-p"]);
		var list = sequences.nucleotideFasta.sequences;
	
		_.each(list, function(seq)  {

	
			//loop through each sequence in the alignment
			seqDinucComposition = {}
		
			var sequence   = seq.sequence;
			var sequenceId = seq.id;
			//glue.log("INFO", "ID result was:", sequenceId);

			//loop through each position in the current sequence
			var lastBase;
			for (var i=0; i < seq.sequence.length; i++) {

				var base = seq.sequence.charAt(i);
			
				if (lastBase) {
			
					var dinuc = lastBase += base;
			
					if (seqDinucComposition[dinuc]) {				
						seqDinucComposition[dinuc] += 1;		
					}
					else {
						seqDinucComposition[dinuc] = 1;
					}				
					lastBase = undefined;
				
				}
				else {
			
				  lastBase = base;
				}		
			}
		
			//glue.log("INFO", "Dinucleotide composition result was:", seqDinucComposition);
			dinucCompositionResults[sequenceId] = seqDinucComposition;
			lengthResults[sequenceId] = seq.sequence.length;
			
		});

 
		_.each(_.keys(dinucCompositionResults), function(sequenceID) {
	
			var seqResults = dinucCompositionResults[sequenceID];
			//glue.log("INFO", "RESULTS FOR '"+sequenceID+"'");

		    var dinucCountResults = {};
		    var dinucFreqResults = {};
            var length = lengthResults[sequenceID];
            
			_.each(_.keys(seqResults), function(dinucleotide) {

				var count = seqResults[dinucleotide];
					
				var ratio = (count / length) * 100;
				var aaFormatedRatio = ratio.toFixed(2);
				//glue.log("INFO", "Got ratio '"+aaFormatedRatio+"' for dinucleotide character'"+dinucleotide+"'");
				dinucCountResults[dinucleotide] = count;
				dinucFreqResults[dinucleotide] = aaFormatedRatio;
							
		
			});

			// add results to array to be returned to GLUE
			outputArray.push({
		
				referenceName: sequenceID,
				seqLength: length,
				AA: dinucFreqResults["AA"],
				AT: dinucFreqResults["AT"],
				AC: dinucFreqResults["AC"],
				AG: dinucFreqResults["AG"],
				TA: dinucFreqResults["TA"],
				TT: dinucFreqResults["TT"],
				TC: dinucFreqResults["TC"],
				TG: dinucFreqResults["TG"],
				CA: dinucFreqResults["CA"],
				CT: dinucFreqResults["CT"],
				CC: dinucFreqResults["CC"],
				CG: dinucFreqResults["CG"],
				GA: dinucFreqResults["GA"],
				GT: dinucFreqResults["GT"],
				GC: dinucFreqResults["GC"],
				GG: dinucFreqResults["GG"]
			   
			});

		
		});


	});

	 
	return outputArray;
}


