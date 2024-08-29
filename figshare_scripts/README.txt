## README 

Supporting data for "Systematic identification of cargo-carrying genetic elements reveals new dimensions of eukaryotic diversity".

# benchmarking/

This directory contains the database of 42 manually annotated elements used for benchmarking the starfish workflow, as well as all other materials required to reproduce the benchmarking results

# starfishScripts.R

This file contains all scripts for generating graphs used in the Figures

# Figure2/

This directory contains data necessary to reproduce Figure 2

# Figure3/

This directory contains data necessary to reproduce Figure 3

# Figure4/

This directory contains data necessary to reproduce Figure 4

# Figure5/

This directory contains data necessary to reproduce Figure 5

# phylo/

This directory contains all sequence alignments and Newick tree files

	/all_superfam_captains.kpicsg.clipkit	trimmed alignment of all predicted Starfish tyrosine recombinase captains
	/all_superfam_captains.kpicsg.clipkit.treefile	a maximum likelihood tree produced from all_superfam_captains.kpicsg.clipkit
	/funTyr50_cap25_crp3.kpicg90.clipkit	trimmed alignment of all mmseqs2 easy-clust representative Starship tyrosine recombinases plus 25 reference Starship tyrosine recombinases plus 3 CryptonF sequences
	/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit	a filtered alignment consisting of the first 512 columns from funTyr50_cap25_crp3.kpicg90.clipkit with all sequences missing 5/6 residues at the predicted active sites removed
	/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile	a maximum likelihood tree produced from funTyr50_cap25_crp3_p1-512_activeFilt.clipkit
	
# sequences/

This directory contains sequences of interest

	/mycodb.final.starships.fna	nucleotide sequences of all predicted Starship elements
	/tyr10771.faa	predicted amino acid sequences of all Starship tyrosine recombinases
	/tyr8369.ids	the sequence IDs of the size-filtered set of 8369 Starship tyrosine recombinases
	

