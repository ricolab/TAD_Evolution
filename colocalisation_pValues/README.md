# TAD_Evolution/colocalisation_pValues
Repository containing code written by Marco Trevisan-Herraz in 2020-2021 to analyse the evidence of the colocalisation of gene ages in pre-defined genomic regions (such as TADs). Code tested in R 3.4.4

Steps to process the data:

1) Copy into the same directory:

	# files with code
	colocalisations_code.R
	colocalisations_functions.R
	
	# files with gene ages
	MiceAges_MTH.tsv
	AgesWithChrX.tsv
	
	# files with genomic regions
	mESCTADs.tsv
	hESCTADs.tsv
	GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_sorted_noNested.txt
	TADs_nCD4_mean_merged.bed
	TADs_Mon_mean_merged.bed
	TADs_Neu_mean_merged.bed

2) Set that folder as working directory using setwd("/myFolder/here/")
3) Run colocalisations_code.R





Note 1: if you want to use your own genomic regions or your own gene ages, you will neet to change those filenames when calling runTADEvolution (the only function within colocalisations_code.R).
Note 2: the randomisation is set using the runif function. In order to get random but reproducible (deterministic) data, the default seed 0 is used in function randomiseAllAges2. If you want to change this, change the seed by calling the function runTADEvolution with an alternative seed using the rndSeed wariable.