# setwd("/pasteHereYourDirectory/")
source("colocalisations_functions.R")
library("biomaRt")

runTADEvolution = function(species,
							analysisType = "tad",
							windowSize = 0,
							chromosomeTable = data.frame(),
							TADFile = "",
							ageFile,
							ageList,
							replacementTable,
							analysisName,
							totRandomisations = 15,
							removePreChr = TRUE,
							rndSeed = 0,
							addRandomisation = TRUE,
							addObserved = TRUE,
							includeLast = TRUE) {
	
	ensemblDataset = "hsapiens_gene_ensembl"
	if(species == "Mus_musculus") ensemblDataset = "mmusculus_gene_ensembl"
	
	analysisType = tolower(analysisType)
	
	if (analysisType == "tad") {
		if (!file.exists(TADFile)) {
			print("Please provide a valid filename when using analysisType == 'TAD'.")
			stop("runTADEvolution stopped.")
			}
		myRegions = read.table(TADFile, header = TRUE)
		}
	if (analysisType == "windows") {
		if (windowSize < 1000) {
			print("Please, include a windowSize >= 1000 when using analysisType == 'TAD'.")
			stop("runTADEvolution stopped.")
			}
		myRegions = generateRegularTADList(chromosomeTable, size = windowSize, includeLast = includeLast)
		}

	if (analysisType != "tad" & analysisType != "windows") {
		print("The parameter analysisType must be either 'Windows' or 'TAD'.")
		stop("runTADEvolution stopped.")
		}

	if (removePreChr) myRegions$chr = substr(myRegions$chr, 4, length(myRegions$chr) - 3)

	ensemblData = useMart(host='feb2014.archive.ensembl.org',
							biomart='ENSEMBL_MART_ENSEMBL',
							dataset=ensemblDataset)
	allBioMart = getBM(mart = ensemblData, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))


	geneAgeList = as.data.frame(read.table(ageFile, header = TRUE, sep = "\t"))
	if (species == "Homo_sapiens") {
		colnames(geneAgeList) = c("ensembl_gene_id", colnames(geneAgeList)[2:6])
		geneAgeList = geneAgeList[,1:6]
		}
	if (species == "Mus_musculus") {
		geneAgeList = as.data.frame(cbind(ensembl_gene_id = as.character(geneAgeList$GeneID), HUGO_symbol = NA, seqnames = NA, start = NA, end = NA, GeneAge = as.character(geneAgeList$GeneAge)))
		}
	
	print(head(myRegions))
	newTADList = getTADGenesTable(allBioMart, myRegions)
	geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeList)
	geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")


	geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]
	geneTADAgeTableAgeOnlyNoDups_YoungMerged = replaceAges(geneTADAgeTableAgeOnlyNoDups, 11, replacementTable)
	

	####### saving output
	ageNumbers_real = ngeneInAge(ageList, geneTADAgeTableAgeOnlyNoDups_YoungMerged)
	write.table(ageNumbers_real, paste0("ageNumbers_", analysisName, "_real.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

	ageNumbers_NH = ngeneInAge(ageList, makeAllOneHugeTAD(geneTADAgeTableAgeOnlyNoDups_YoungMerged))
	write.table(ageNumbers_NH, paste0("ageNumbers_", analysisName, "_NH.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
	
	if(addObserved) {
		myFinalList_real = getAllPairData(geneTADAgeTableAgeOnlyNoDups_YoungMerged, ageList, totRandomisations = totRandomisations)
		MLMatrix_real = fromFinalListToSgnTable(myFinalList_real)
		write.table(MLMatrix_real, paste0("minusLogPValueSgn_table_", analysisName, "_real.tsv"), sep = "\t")

		MLMatrix_NH = fromFinalListToSgnTable(myFinalList_real, getHugeTADdataInstead = TRUE)
		write.table(MLMatrix_NH, paste0("minusLogPValueSgn_table_", analysisName, "_NH.tsv"), sep = "\t")
		}
	
	# rndAge part, scrambling age labels
	if(addRandomisation) {
		geneTADAgeTableAgeOnlyNoDups_YoungMerged_rndAge = randomiseAllAges(geneTADAgeTableAgeOnlyNoDups_YoungMerged, "gene_age", seed = rndSeed)
		myFinalList_rndAge = getAllPairData(geneTADAgeTableAgeOnlyNoDups_YoungMerged_rndAge, ageList, totRandomisations = totRandomisations)
		MLMatrix_rndAge = fromFinalListToSgnTable(myFinalList_rndAge)
		write.table(MLMatrix_rndAge, paste0("minusLogPValueSgn_table_", analysisName, "_rndAge.tsv"), sep = "\t")
		}
	
	}


# mESC merging GliresRodentia
runTADEvolution(
			species = "Mus_musculus",
			TADFile = "mESCTADs.tsv",
			ageFile = "MiceAges_MTH.tsv", # from Caelinn's email of 16mar-2020, file "MiceAges.txt", changing only Fungi-Metazoa_group to FungiMetazoa
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC"
			)
			
# hESC merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "hESCTADs.tsv",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC"
			)
			
# Bcell merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "GSE63525_GM12878_primary+replicate_Arrowhead_domainlist_sorted_noNested.txt",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hBcell",
			removePreChr = FALSE
			)

# nCD4 merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "TADs_nCD4_mean_merged.bed",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "nCD4",
			removePreChr = FALSE
			)

# monocytes merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "TADs_Mon_mean_merged.bed",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "mono",
			removePreChr = FALSE
			)

# neutrophils merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "TADs_Neu_mean_merged.bed",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "neut",
			removePreChr = FALSE
			)

#
#
#
# analyses with regular windows
			
# hESC regular 5 kb
runTADEvolution(
			species = "Homo_sapiens",
			analysisType = "windows",
			windowSize = 5000,
			chromosomeTable = chromosomes_hg19,
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC_5kb"
			)

# hESC regular 50 kb
runTADEvolution(
			species = "Homo_sapiens",
			analysisType = "windows",
			windowSize = 50000,
			chromosomeTable = chromosomes_hg19,
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC_50kb"
			)

# hESC regular 500 kb
runTADEvolution(
			species = "Homo_sapiens",
			analysisType = "windows",
			windowSize = 500000,
			chromosomeTable = chromosomes_hg19,
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC_500kb"
			)
			
# hESC regular 5 Mb
runTADEvolution(
			species = "Homo_sapiens",
			analysisType = "windows",
			windowSize = 5000000,
			chromosomeTable = chromosomes_hg19,
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC_5000kb"
			)
			
# mESC regular 5 kb
runTADEvolution(
			species = "Mus_musculus",
			analysisType = "windows",
			windowSize = 5000,
			chromosomeTable = chromosomes_mm10,
			ageFile = "MiceAges_MTH.tsv",
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC_5kb"
			)
			
# mESC regular 50 kb
runTADEvolution(
			species = "Mus_musculus",
			analysisType = "windows",
			windowSize = 50000,
			chromosomeTable = chromosomes_mm10,
			ageFile = "MiceAges_MTH.tsv",
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC_50kb"
			)

# mESC regular 500 kb
runTADEvolution(
			species = "Mus_musculus",
			analysisType = "windows",
			windowSize = 500000,
			chromosomeTable = chromosomes_mm10,
			ageFile = "MiceAges_MTH.tsv",
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC_500kb"
			)
			
# mESC regular 5 Mb
runTADEvolution(
			species = "Mus_musculus",
			analysisType = "windows",
			windowSize = 5000000,
			chromosomeTable = chromosomes_mm10,
			ageFile = "MiceAges_MTH.tsv",
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC_5000kb"
			)
			
#
#
# some more randomisations for mESC and hESC

runTADEvolution(
			species = "Mus_musculus",
			TADFile = "mESCTADs.tsv",
			ageFile = "MiceAges_MTH.tsv", # from Caelinn's email of 16mar-2020, file "MiceAges.txt", changing only Fungi-Metazoa_group to FungiMetazoa
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC_seed1",
			addObserved = FALSE,
			rndSeed = 1
			)

# hESC merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "hESCTADs.tsv",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC_seed1",
			addObserved = FALSE,
			rndSeed = 1
			)

runTADEvolution(
			species = "Mus_musculus",
			TADFile = "mESCTADs.tsv",
			ageFile = "MiceAges_MTH.tsv", # from Caelinn's email of 16mar-2020, file "MiceAges.txt", changing only Fungi-Metazoa_group to FungiMetazoa
			ageList = ageListMouse,
			replacementTable = replacementTableMouse_GliresRodentia,
			analysisName = "mESC_seed2",
			addObserved = FALSE,
			rndSeed = 2
			)

# hESC merging Primates
runTADEvolution(
			species = "Homo_sapiens",
			TADFile = "hESCTADs.tsv",
			ageFile = "AgesWithChrX.tsv", # converted from Caelinn's file in /data
			ageList = ageListHuman,
			replacementTable = replacementTableHuman,
			analysisName = "hESC_seed2",
			addObserved = FALSE,
			rndSeed = 2
			)