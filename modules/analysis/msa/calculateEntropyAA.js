function calculateEntropyAA(refsequence, alignment, source, feature) {

	// set possible amino acids
	// check each aa is valid?
	//var expectedAminoAcids = ['A', 'I', 'L', 'M', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P', '*'];

	var alignmentString = "alignment/" + alignment;
	var sourceName = "sequence.source.name = '" + source + "'"
	// get frequencies for chosen alignment and store in object
	values_obj = {};
	glue.inMode(alignmentString, function(){
	
		var aaFreqs = glue.command(["amino-acid", "frequency", "-w", sourceName, "-r", refsequence, "-f", feature]);
		var info = aaFreqs["alignmentAminoAcidFrequencyResult"]["row"];
		_.each(info, function(count){
	
			var codon_pos = count["value"][1];
			var aminoAcid = count["value"][2];
			var numMembers = count["value"][3];
				
			if (values_obj[codon_pos] == null){
				values_obj[codon_pos] = [{aminoAcid : aminoAcid, count : numMembers }];
			} else {
				var val = {aminoAcid : aminoAcid, count : numMembers };
				values_obj[codon_pos].push(val);
			}
		});	
	});

	//calculate entropy for each codon position
	var position = 1;
	var entropy_array = [];
	_.each(values_obj, function(pos){

		var entropy = 0;
		var counts = [];
		var total = 0;
	
		//extract counts and total at current position
		_.each(pos, function(aa){
			var num = aa['count'];
			total += num;
			counts.push(num);		
		});
		// generate values for each aa for entropy calculation
		var entropy_values = [];
		_.each(counts, function(val){
			ev = -(val /total) * Math.log(val / total);
			entropy_values.push(ev);
		});
		// calculate entropy
		_.each(entropy_values, function(ev){
			entropy += ev;
		});
		// add results to array to be returned to GLUE
		entropy_array.push({
			codon: position,
			entropy: entropy.toFixed(3)
		});
		position++;
	});

	//output result array
	return entropy_array;

}
