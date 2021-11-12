# added on 16nov-2020 to include rndAge
randomiseAllAges = function(myTable, ageColumn = "gene_age", chrColumn = "gene_chr", seed = -1, randomisationType = "full") {

	# randomisationType can be
	#	"full" (default, scrambling the whole table, keeping the age frequencies)
	#	"chr" (scrambling within each chromosome, keeping  the age frequencies within each chromosome)
	# Warning! myTable should be sorted in advance by {gene_chr, gene_start, gene_end} to make it work in chr mode.
	
	# use seed != -1 if you don't want a deterministic random distribution
	if (seed != -1) set.seed(seed)
	newTable = myTable
	myAges = as.character(myTable[, ageColumn])
	myChrs = as.character(myTable[, chrColumn])
	scramblingTable = as.data.frame(cbind(myChrs, myAges, runif(length(myAges))))
	colnames(scramblingTable) = c(chrColumn, ageColumn, "rnd")
	
	if (randomisationType == "full") scramblingTable = scramblingTable[order(scramblingTable$rnd),]
	if (randomisationType == "chr") scramblingTable = scramblingTable[order(scramblingTable[, chrColumn], scramblingTable$rnd),]
	
	newTable[, ageColumn] = scramblingTable[, ageColumn]
	
	return(newTable)
	}
	
saveRandomFiles = function(ageFile, amount = 30, ageColumn = "gene_age", chrColumn = "gene_chr", startColumn = "start", endColumn = "end", idColumn = "ensembl_gene_id") {
	totRandomisations = amount
	geneAgeList = as.data.frame(read.table("MiceAges_MTH_withCoord.tsv", header = TRUE, sep = "\t"))
	geneAgeList = geneAgeList[c(idColumn, chrColumn, startColumn, endColumn, ageColumn)]
	geneAgeList = geneAgeList[order(geneAgeList[, chrColumn], geneAgeList[, startColumn], geneAgeList[, endColumn]),]
	colnames(geneAgeList) = c(idColumn, chrColumn, startColumn, endColumn, ageColumn)
		
	for (n in 1:totRandomisations) {
		randTable = randomiseAllAges(geneAgeList, ageColumn = ageColumn, chrColumn = chrColumn, seed = n, randomisationType = "chr")
		newFileName = paste0(ageFile, "_rnd", sprintf("%02d",n), "tsv")
		colnames(randTable) = c(idColumn, chrColumn, startColumn, endColumn, ageColumn)
		write.table(randTable, newFileName, row.names = FALSE, quote = FALSE)
		}
	}

# Primates = Simiiformes + Catarrhini + Hominoidea + Hominidae + HomoPanGorilla + HomoSapiens
# Eutheria = Eutheria + Glires
# Rodentia = Rodentia + Sciurognathi + Murinae + Mus.musculus

replaceAges = function(myTable, myColumnNumber, replacementTable) {
    
    resultingTable = myTable
	resultingTable[,myColumnNumber] = as.character(resultingTable[,myColumnNumber])
	replacementTable[,1] = as.character(replacementTable[,1])
	replacementTable[,2] = as.character(replacementTable[,2])
	
    for (i in 1:nrow(replacementTable)) {
		resultingTable[resultingTable[,myColumnNumber] == replacementTable[i, "from"], colnames(resultingTable)[myColumnNumber]] = replacementTable[i, "to"]
        }
		
    return(resultingTable)
    
    }
    
# this function is similar to fillTableWithGeneTypes
# input myTADs have four columns:
#     chr --> as 1, 2, 3, ... X
#     start
#     end
#     percentOverlap
# input myGenes have six columns:
#     ensembl_gene_id
#     wikigene_name
#     chromosome_name --> as 1, 2, 3, ... X
#     start_position
#     end_position
#     gene_biotype
#

library(ggplot2)

filterTableWithList = function(myGeneTADAgeTable, myList, filtering = "inclusion") {
	
	if(filtering != "inclusion" & filtering != "exclusion") {
		print("Error, filtering can only be 'inclusion' or 'exclusion'.")
		return(NULL)
		}
		
	if (filtering == "inclusion") newGeneTADAgeTable = subset(myGeneTADAgeTable, ensembl_gene_id %in% unlist(myList))
	if (filtering == "exclusion") newGeneTADAgeTable = subset(myGeneTADAgeTable, !(ensembl_gene_id %in% unlist(myList)))
	
	return(newGeneTADAgeTable)
	}

fromCSVAgesToTSVAges = function(inputFile, outputFile) {
	# tsv file "/mnt/sexreg/temporal/geneage/pgen.1007902.s021.tsv"
	# csv file "/mnt/sexreg/temporal/geneage/Ages.csv"
	csvTable = read.csv(file = inputFile, header = TRUE)
	tsvTable = cbind(as.character(csvTable$X), rep("-", nrow(csvTable)), as.character(csvTable$Chromosome), as.numeric(as.character(csvTable$Start)), as.numeric(as.character(csvTable$End)), as.character(csvTable$GeneAge), as.character(csvTable$BirthType))
	colnames(tsvTable) = c("X", "HUGO_symbol", "seqnames", "start", "end", "GeneAge", "BirthType")
	
	write.table(tsvTable, file = outputFile, sep = "\t", row.names = FALSE)
	print("Done.")
	}

getTADGenesTable = function(myGenes, myTADs) {
	feedbackLoop = 1000
	# all.equal(newTADList_old, newTADList) gives true
	# without doing this stack thing, it took 6:15
	# using stacks of 3000, it takes 1:17 to do the same
	# using stacks of 100, it takes 1:13
	stackSize = 100
	TADGenesTable = data.frame()
	TADGenesTableStack = data.frame()
	geneDistribution = c()
	for (i in 1:nrow(myGenes)) {
		if (i %% feedbackLoop == 0) print(paste0("Checking gene #", i, "/", nrow(myGenes)))
		row = myGenes[i,]
		mySubset = subset(myTADs, myTADs$end > row$start_position & myTADs$start < row$end_position & myTADs$chr == row$chromosome_name)
		geneDistribution = c(geneDistribution, nrow(mySubset))
		if (nrow(mySubset) == 1) {
			newRow = cbind(row, mySubset[, 1:3], paste(mySubset[1,1], mySubset[1,2], mySubset[1,3], sep = "_"))
			TADGenesTableStack = rbind(TADGenesTableStack, newRow)
			}
			
		if (i %% stackSize == 0 | i == nrow(myGenes)) {
			TADGenesTable = rbind(TADGenesTable, TADGenesTableStack)
			TADGenesTableStack = data.frame()
			}
		}
		
		for (j in 0:max(geneDistribution)) print(paste0("There are ", sum(geneDistribution == j), " genes present in exactly ", j, " TADs (", round(sum(geneDistribution == j) / length(geneDistribution) * 100, 2), "%)."))
		
		colnames(TADGenesTable) = c("ensembl_gene_id", "wikigene_name", "gene_chr", "gene_start", "gene_end", "gene_biotype", "tad_chr", "tad_start", "tad_end", "tad_id")
		
		return(TADGenesTable)
	}
	
# this function uses two input files:
# 1) the output of getTADGenesTable,
# 2) a file including ensembl gene ids with ages, 
# which is like this one downloadable from https://ndownloader.figshare.com/files/14163824
# in turn, published here https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007902#sec025
# the columns for this file are (after removing the last columns with geneAgeListHuman = geneAgeListHuman[,1:6]):
# (unnamed) ensembl_gene_id
# HUGO symbol
# seqnames --> the chromosome_name
# start
# end
# GeneAge, which in turn can be any of:
#      HomoSapiens
#      Euteleostomi
#      HomoPanGorilla
#      Tetrapoda
#      Eutheria
#      Bilateria
#      Amniota
#      Chordata
#      FungiMetazoa
#      Theria
#      Hominoidea
#      Sarcopterygii
#      Hominidae
#      Simiiformes
#      Catarrhini
#      Mammalia
# from these, Caelinn has merged into "Primate":
#      Simiiformes
#      Catarrhini
#      Hominoidea
#      Hominidae
#      HomoPanGorilla
#      HomoSapiens

mergeGeneTADAge = function(myGeneTADs, myGeneAges) {
	
	colnames(myGeneAges) = c("ensembl_gene_id", colnames(myGeneAges)[2:6])
	gene_age = c()
	for (i in 1:nrow(myGeneTADs)) {
		# if there are repeated ensembl_gene_id elements, this will give problems
		thisAge = as.character(subset(myGeneAges, as.character(ensembl_gene_id) == myGeneTADs[i, 1])[1,"GeneAge"])
		if (is.na(thisAge)) thisAge = "No Age Provided"
		gene_age = c(gene_age, thisAge)
		}

	geneTADAge = cbind(myGeneTADs, gene_age)
	return(geneTADAge)
	}
	
# must be a table with a column called tad_id, such as the result of
# getTADGenesTable or mergeGeneTADAge
ntad = function(myTable) {
	myTADids = as.data.frame(myTable$tad_id)
	myTADids = unique(myTADids)
	
	return(nrow(myTADids))
	}
	
ngene = function(myTable, myTADid) {
	TADrows = unique(subset(myTable, tad_id == myTADid))
	return(nrow(TADrows))
	}
	
getagelist = function(myTable) {
	myAgeList = unique(myTable[, "gene_age"])
	return(myAgeList)
	}
	
getTADsWithAgePerChromosome = function(myTable, myAgeList) {
	thisTable = data.frame()
	for (chr in c(1:22, "X", "Y")) {
		chrName = paste0("chr", chr)
		totTADs = ntad(subset(myTable, tad_chr == chr))
		chrRow = c(chrName, totTADs)
		for (age in as.character(unlist(myAgeList))) {
			totTADsInAge = ntad(subset(myTable, tad_chr == as.character(chr) & gene_age == age))
			chrRow = c(chrRow, totTADsInAge)
			}
		chrRow = t(as.data.frame(chrRow))
		thisTable = rbind(thisTable, chrRow)
		}
		
	colnames(thisTable) = c("chr", "totTADs", as.character(unlist(myAgeList)))
	rownames(thisTable) = NULL
	return(thisTable)
	}

getColocalisationProbability = function(myTable, agePre, agePost, makePlot = FALSE, outputOnlyProbability = TRUE, verbose = TRUE) {
	# first, I get the list of all the genes with that age
	myListOfGenes = unique(subset(myTable, gene_age == agePre)[, "ensembl_gene_id"])
	if (verbose) print(paste0("Total number of genes with age ", agePre, ": ", length(myListOfGenes)))
	probList = data.frame()
	coutner = 0
	
	if (length(myListOfGenes) == 0) {
		probList = rbind(probList, t(as.data.frame(c(NA, NA, NA))))
		}
	else {
		for (thisGene in myListOfGenes) {
			coutner = coutner + 1
			if (coutner %% 200 == 0 & verbose) print(paste0(coutner, "..."))
			relevantTAD = subset(myTable, ensembl_gene_id == thisGene)[, "tad_id"]
			geneList = subset(myTable, tad_id == relevantTAD & ensembl_gene_id != thisGene)
			totGenesInTAD = nrow(geneList)
			totGenesInTADWithAgePost = nrow(subset(myTable, tad_id == relevantTAD & ensembl_gene_id != thisGene & gene_age == agePost))
			if (totGenesInTAD > 0) {
				probList = rbind(probList, t(as.data.frame(c(thisGene, as.character(relevantTAD), as.numeric(totGenesInTADWithAgePost / totGenesInTAD)))))
				}
			}
		}
	if (nrow(probList) == 0) probList = t(as.data.frame(c(NA, NA, NA)))
	colnames(probList) = c("ensembl_gene_id", "tad_id", "probability")
	rownames(probList) = c()
	probList[, "probability"] = as.numeric(as.character((probList[, "probability"])))
	
	if (makePlot) {
		thisPlot = ggplot(probList, aes(x = probability)) + geom_histogram(aes(y = ..density.., colour = "black", fill = "white"), colour = "black", fill = "white", bins = 20) + geom_density(alpha = 0.2, fill = "#40ff40") + xlim(c(0, 1)) + ylim(c(0,6))
		# print(thisPlot)
		}
	overallProbability = mean(probList[, "probability"], na.rm = TRUE)
	if (verbose) print(paste0("Overall probability: ", overallProbability))
	
	if (outputOnlyProbability) return(overallProbability)
	else return(probList)
	}

interTADRemover = function(myTADs, myChromosomes) {
	
	# indices for myChromosomes
	CHROMINDEX = 1
	CENTROMEREINDEX = 2
	CHRENDINDEX = 3
	
	# indices for myTAD
	CHROMINDEX = 1
	STARTINDEX = 2
	ENDINDEX = 3
	TADNAMEINDEX = 4
	
	chromList = myChromosomes[,CHROMINDEX]
	
	newTADs = vector()
	for (chromNum in 1:length(chromList)) {
		thisChrom = chromList[chromNum]
		thisChromNoPrefix = substr(thisChrom, 4, 5)
		
		print(paste0("Working on ", thisChrom, "..."))
		TADsubset = as.data.frame(subset(myTADs, paste0("chr", chr) == thisChrom))
		if (nrow(TADsubset) > 0) {
			
			startingPoint = 0
			for (i in 1:nrow(TADsubset)) {
				
				TADName = paste0("Original tad_id = ", TADsubset$chr[[i]], "_", TADsubset$start[[i]], "_", TADsubset$end[[i]])
				
				if (i < nrow(TADsubset)) endingPoint = floor((TADsubset$end[[i]] + TADsubset$start[[i + 1]]) / 2)
				else endingPoint = as.numeric(as.character(subset(as.data.frame(myChromosomes), chrom == thisChrom)$end))
				
				newRow = c(chr = thisChromNoPrefix, start = startingPoint, end = endingPoint, TAD = TADName, length = endingPoint - startingPoint)
				
				newTADs = rbind(newTADs, newRow)
				
				startingPoint = endingPoint
				
				}
			
			}
		else print(paste0("Warning: no TADs for ", thisChrom, ", removing..."))
		}
		
			
	rownames(newTADs) = NULL
	newTADs = as.data.frame(newTADs)
	newTADs$chr = as.character(newTADs$chr)
	newTADs$start = as.numeric(as.character(newTADs$start))
	newTADs$end = as.numeric(as.character(newTADs$end))
	newTADs$TAD = as.character(as.character(newTADs$TAD))
	newTADs$length = as.numeric(as.character(newTADs$length))
	return(newTADs)
	
	}

fusedTAD = function(myTADs, myChromosomes) {
# forked from version of 12dic-2019
# note: this algorithm removes the interTAD regions
	
	# indices for myChromosomes
	CHROMINDEX = 1
	CENTROMEREINDEX = 2
	CHRENDINDEX = 3
	
	# indices for myTAD
	CHROMINDEX = 1
	STARTINDEX = 2
	ENDINDEX = 3
	TADNAMEINDEX = 4
	
	chromList = myChromosomes[,CHROMINDEX]
	
	newTADs = vector()
	for (chromNum in 1:length(chromList)) {
		thisChrom = chromList[chromNum]
		thisChromNoPrefix = substr(thisChrom, 4, 5)
		
		print(paste0("Working on ", thisChrom, "..."))
		TADsubset = as.data.frame(subset(myTADs, paste0("chr", chr) == thisChrom))
		if (nrow(TADsubset) > 0) {
		
			startingPoint = 0
			previousTAD = paste0("beginning of ", thisChrom)
			for (i in 1:nrow(TADsubset)) {
				# print(paste0(thisChrom, ":", TADsubset$start[[i]], "-", TADsubset$end[[i]]))
				cuttingPoint = floor((TADsubset$start[[i]] + TADsubset$end[[i]]) / 2)
				nextTAD = paste0(as.character(TADsubset[i, CHROMINDEX]), "_", as.character(TADsubset[i, STARTINDEX]), "_", as.character(TADsubset[i, ENDINDEX]))
				previousTADMessage = previousTAD
				nextTADMessage = nextTAD
				if (previousTAD != paste0("beginning of ", thisChrom)) previousTADMessage = paste0("second half of ", previousTADMessage)
				if (nextTAD != paste0("ending of ", thisChrom)) nextTADMessage = paste0("first half of ", nextTADMessage)
				newRow = c(chr = thisChromNoPrefix, start = startingPoint, end = cuttingPoint, TAD = paste0(previousTADMessage, " & ", nextTADMessage), length = cuttingPoint - startingPoint)
				
				newTADs = rbind(newTADs, newRow)
				
				previousTAD = nextTAD
				startingPoint = cuttingPoint
				}
				
			# now adding last half TAD
				cuttingPoint = as.numeric(as.character(subset(as.data.frame(myChromosomes), chrom == thisChrom)$end))
				nextTAD = paste0("ending of ", thisChrom)
				previousTADMessage = previousTAD
				nextTADMessage = nextTAD
				if (previousTAD != paste0("beginning of ", thisChrom)) previousTADMessage = paste0("second half of ", previousTADMessage)
				if (nextTAD != paste0("ending of ", thisChrom)) nextTADMessage = paste0("first half of ", nextTADMessage)
				newRow = c(chr = thisChromNoPrefix, start = startingPoint, end = cuttingPoint, TAD = paste0(previousTADMessage, " & ", nextTADMessage), length = cuttingPoint - startingPoint)
				
				newTADs = rbind(newTADs, newRow)
		}
		else print(paste0("Warning: no TADs for ", thisChrom, ", removing..."))
		
		}
		
	rownames(newTADs) = NULL
	newTADs = as.data.frame(newTADs)
	newTADs$chr = as.character(newTADs$chr)
	newTADs$start = as.numeric(as.character(newTADs$start))
	newTADs$end = as.numeric(as.character(newTADs$end))
	newTADs$TAD = as.character(as.character(newTADs$TAD))
	newTADs$length = as.numeric(as.character(newTADs$length))
	return(newTADs)
	
	}
# fff = fusedTAD(TADList, chromosomes_hg19)

makeAllOneHugeTAD = function(myTable, hugeTADName = "hugeTAD") {
	newTable = myTable
	newTable$tad_id = hugeTADName
	
	return(newTable)
	}
	

getPairData = function(TADTable, age1 = "", age2 = "", totRandomisations = 30, verbose = FALSE, showGraph = FALSE, randomisationType = "full") {
	maxSgnMinusLofPValue = 30
	
	if(nchar(age1) == 0 | nchar(age2) == 0) stop("Please don't leage age1 or age2 empty.")
	
	print("")
	print("********************")
	print(paste0("********** age1: ", age1, ", age2: ", age2))
	print("Calculating real probability.")
	realProbability = getColocalisationProbability(TADTable, age1, age2, verbose = verbose)
	
	if (is.na(realProbability) | is.nan(realProbability)) {
		print(paste0("There was a problem calculating the probability, probably there are no genes in this table with age ", age1))
		finalList = c(age1, age2, NA, NA, NA, NA, NA, NA, rep(NA, totRandomisations))
		}
	else {
		print ("Calculating huge TAD probability.")
		hugeTADProbability = getColocalisationProbability(makeAllOneHugeTAD(TADTable), age1, age2, verbose = verbose)
		randomProbabilityList = c()
		for (i in 1:totRandomisations) {
			print (paste0("Calculating randomisation #", i, "/", totRandomisations, "."))
			randomProbability = getColocalisationProbability(randomiseAllAges(TADTable, "gene_age", randomisationType = randomisationType), age1, age2, verbose = verbose)
			randomProbabilityList = c(randomProbabilityList, randomProbability)
			}
		randomProbabilityList = sort(randomProbabilityList, na.last = TRUE)
		meanRandomisations = mean(randomProbabilityList, na.rm = TRUE)
		sdevRandomisations = sd(randomProbabilityList, na.rm = TRUE)
		
		dfrandomProbabilityList = as.data.frame(randomProbabilityList)
		colnames(dfrandomProbabilityList) = "prob"
		
		if (showGraph) {
			arrow = arrow(angle = 15, type = "closed")
			myPlot = ggplot(dfrandomProbabilityList, aes(x = 1, y = prob)) + geom_boxplot() + geom_jitter() + ylim(0,1) + xlim(0,10) + annotate("segment", x = 0.5, xend = 1, y = realProbability, yend = realProbability, arrow = arrow, colour = "red") + annotate("segment", x = 0.5, xend = 1, y = hugeTADProbability, yend = hugeTADProbability, arrow = arrow, colour = "black")
			print(myPlot)
			}
		
		if (realProbability <= meanRandomisations) realPValue = 2 * pnorm(realProbability, mean = meanRandomisations, sd = sdevRandomisations)
		else realPValue = 2 * (1 - pnorm(realProbability, mean = meanRandomisations, sd = sdevRandomisations))
		if (hugeTADProbability <= meanRandomisations) hugeTADPValue = 2 * pnorm(hugeTADProbability, mean = meanRandomisations, sd = sdevRandomisations)
		else hugeTADPValue = 2 * (1 - pnorm(hugeTADProbability, mean = meanRandomisations, sd = sdevRandomisations))
		
		if (realPValue == 0) sgnMinusLogPValueReal = maxSgnMinusLofPValue * sign(realProbability - meanRandomisations)
		else sgnMinusLogPValueReal = min(-log(realPValue), maxSgnMinusLofPValue) * sign(realProbability - meanRandomisations)
		
		if (hugeTADPValue == 0) sgnMinusLogPValueHugeTAD = maxSgnMinusLofPValue * sign(realProbability - meanRandomisations)
		else sgnMinusLogPValueHugeTAD = min(-log(hugeTADPValue), maxSgnMinusLofPValue) * sign(realProbability - meanRandomisations)
		
		print ("")
		print(paste0("Real probability: ", realProbability, ", p-value: ", realPValue))
		print(paste0("Huge TAD probability: ", hugeTADProbability, ", p-value: ", hugeTADPValue))
		print("Random probability distribution:")
		if(length(randomProbabilityList) != totRandomisations) print(paste0("WARNING: missing randomiosations: randomProbabilityList = ", length(randomProbabilityList), ", totRandomisations = ", totRandomisations))
		for (i in 1:length(randomProbabilityList)) print(paste0(i, ": ", randomProbabilityList[i]))
		print("")
		print(paste0("Mean randomisations: ", meanRandomisations))
		print(paste0("Sigma randomisations: ", sdevRandomisations))
		
		finalList = c(age1, age2, realProbability, realPValue, sgnMinusLogPValueReal, hugeTADProbability, hugeTADPValue, sgnMinusLogPValueHugeTAD, randomProbabilityList)
		}
	dffinalList = as.data.frame(finalList)
	rownames(dffinalList) = c("age1", "age2", "probReal", "pValueReal", "minusLogPvalueSgnReal", "probHugeTAD", "pValueHugeTaD", "minusLogPvalueSgnHugeTAD", paste0("rnd", c(1:totRandomisations)))
	colnames(dffinalList) = paste0(age1, "/", age2)
	
	return(dffinalList)
	}
	
showBoxPlot = function(dataTable) {
	POSREALPROBABILITY = 3
	POSHUGETADPROBABILITY = 6
	POSPVALUEHUGETAD = 7
	POSLOGPVALUEHUGETAD = 8
	
	arrow = arrow(angle = 15, type = "closed")
	randomData = data.frame()
	for (i in 1:ncol(dataTable)) {
		firstColumn = rep(i, nrow(dataTable) - POSLOGPVALUEHUGETAD)
		secondColumn = as.numeric(as.character(dataTable[(POSLOGPVALUEHUGETAD + 1):nrow(dataTable), i]))
		subTable = as.data.frame(cbind(firstColumn, secondColumn))
		randomData = rbind(randomData, subTable)
		}
	colnames(randomData) = c("ageCombination", "prob")
	myPlot = ggplot(randomData, aes(x = ageCombination, y = prob, group = ageCombination)) + geom_boxplot() + ylim(0,1) + xlim(0, ncol(dataTable) + 1) + geom_jitter()

	for (i in 1:ncol(dataTable)) {
		realProbability = as.numeric(as.character(dataTable[POSREALPROBABILITY, i]))
		myPlot = myPlot + annotate("segment", x = i - 0.5, xend = i, y = realProbability, yend = realProbability, arrow = arrow, colour = "red")
		}
		
	for (i in 1:ncol(dataTable)) {
		hugeTADProbability = as.numeric(as.character(dataTable[POSHUGETADPROBABILITY, i]))
		myPlot = myPlot + annotate("segment", x = i - 0.5, xend = i, y = hugeTADProbability, yend = hugeTADProbability, arrow = arrow, colour = "black")
		}
		
	print(myPlot)
	
	}

getManyPairData = function(TADTable, specificAge, ageList, totRandomisations = 30, verbose = FALSE, showGraph = FALSE, randomisationType = "full") {
	finalTable = data.frame()
	for(i in 1:nrow(ageList)) {
		thisAge = as.character(ageList[i,])
		colAge = getPairData(TADTable, specificAge, as.character(thisAge), totRandomisations = totRandomisations, verbose = verbose, randomisationType = randomisationType)
		if(nrow(finalTable) == 0) finalTable = colAge
		else finalTable = cbind(finalTable, colAge)
		}
		
	if(showGraph) showBoxPlot(finalTable)
		
	return(finalTable)
	
	}
	
getAllPairData = function(TADTable, ageList, totRandomisations = 30, verbose = FALSE, showGraph = FALSE, randomisationType = "full") {
	
	resultingList = list()
	
	for(i in 1:nrow(ageList)) {
		thisAge = as.character(ageList[i,])
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print(paste0("++++++++++++++++++++ Starting with interactions for: ", thisAge, " ++++++++++++++++++++"))
		dfAge = getManyPairData(TADTable, thisAge, ageList, totRandomisations = totRandomisations, verbose = verbose, showGraph = showGraph, randomisationType = randomisationType)
		resultingList[[i]] = dfAge
		}
		
	return(resultingList)
	}
	
fromFinalListToSgnTable = function(finalList, getHugeTADdataInstead = FALSE) {
	POSAGE1 = 1
	POSAGE2 = 2
	POSMINUSLOGREAL = 5
	POSMINUSLOGHUGETAD = 8
	
	tableName = "minusLogReal"
	if (getHugeTADdataInstead) tableName = "minusLogHugeTAD"
	
	minusLogMatrix = data.frame()
	for (i in 1:length(finalList)) {
		if (getHugeTADdataInstead) minusLogRow = as.data.frame(finalList[[i]][POSMINUSLOGHUGETAD,])
		else minusLogRow = as.data.frame(finalList[[i]][POSMINUSLOGREAL,])
		colnames(minusLogRow) = as.character(unlist(finalList[[i]][POSAGE2,]))
		rownames(minusLogRow) = as.character(unlist(finalList[[i]][POSAGE1,1]))

		minusLogMatrix = rbind(minusLogMatrix, minusLogRow)
		}
	
	return(minusLogMatrix)
	
	}
	
fromFinalListToProbDiffTable = function(finalList) {
	POSAGE1 = 1
	POSAGE2 = 2
	POSPROB = 3
	POSPVALUE = 4
	POSMINUSLOGREAL = 5
	POSPROBNH = 6
	POSPVALUENH = 7
	POSMINUSLOGHUGETAD = 8
	
	probDiffMatrix = data.frame()
	for (i in 1:length(finalList)) {
		probsReal = as.numeric(as.character(unlist(finalList[[i]][POSPROB,])))
		probsNH = as.numeric(as.character(unlist(finalList[[i]][POSPROBNH,])))
		
		probDiffRow = t(as.data.frame((probsReal - probsNH) / probsNH))
		colnames(probDiffRow) = as.character(unlist(finalList[[i]][POSAGE2,]))
		rownames(probDiffRow) = as.character(unlist(finalList[[i]][POSAGE1,1]))

		probDiffMatrix = rbind(probDiffMatrix, probDiffRow)
		}
	
	return(probDiffMatrix)
	
	}
# probDiffMatrix_14abr_mESC_iTR = fromFinalListToProbDiffTable(myFinalList_14abr_mESC_iTR)

ngeneInAge = function(myAgeList, myGeneTADAgeTable) {
	totGenes = nrow(myGeneTADAgeTable)
	outputDF = data.frame()
	for (agestr in myAgeList$age) {
		ngenes = nrow(subset(myGeneTADAgeTable, as.character(gene_age) == agestr))
		cat(paste0(agestr, ",", ngenes, ",", ngenes / totGenes, "\n"))
		thisRow = t(data.frame(c(age = agestr, ngenes = ngenes, ratio = ngenes / totGenes)))
		rownames(thisRow) = NULL
		outputDF = rbind(outputDF, thisRow)
		}
		
	return(outputDF)
	}

generateRegularTADList = function(chromosomes, size = 500000, includeLast = TRUE) {
	
	finalList = data.frame()
	for (i in 1:nrow(chromosomes)) {
		chrom = chromosomes[i, "chrom"]
		chromEnd = as.numeric(chromosomes[i, "end"])
		chromRows = ceiling(chromEnd / size)
		chromTotal = data.frame(chr = rep(chrom, chromRows),
								start = (0:(chromRows - 1)) * size,
								end = (1:chromRows) * size - 1)
		if (chromTotal[chromRows,"end"] > chromEnd) {
			if (includeLast) chromTotal[chromRows,"end"] = chromEnd
			else chromTotal = head(chromTotal, -1)
			}
		finalList = rbind(finalList, chromTotal)
		}
	
	rownames(finalList) = NULL
	finalList$chr = as.character(unlist(finalList$chr))
	finalList$start = format(as.numeric(as.character(unlist(finalList$start))), scientific = FALSE)
	finalList$end = format(as.numeric(as.character(unlist(finalList$end))), scientific = FALSE)
	
	return(finalList)
	
	}
	

################################################ functions for essentiality

getGeneEssentialityTable = function(geneList, essentialGenes) {
	
	# ensemble_gene_id = igual
	# HUGO_symbol -> "-"
	# seqnames = chromosome_name
	# start -> start_position
	# end -> end_position
	# GeneAge -> essentiality
	
	geneEssList = cbind(ensembl_gene_id = geneList[,"ensembl_gene_id"], HUGO_symbol = "-", seqnames = geneList[,"chromosome_name"], start = geneList[,"start_position"], end = geneList[,"end_position"], GeneAge = "nonEssential")
	geneEssList = as.data.frame(geneEssList, stringsAsFactors = FALSE)
	for (i in 1:nrow(essentialGenes)) {
		if (i %% 1000 == 0) print(paste0("Checked ", i, "/", nrow(essentialGenes), " rows."))
		
		thisGene = essentialGenes[i,1]
		geneEssList[geneEssList$ensembl_gene_id == thisGene, "GeneAge"] = "essential"
		}

	rownames(geneEssList) = NULL
	return(geneEssList)
	
	}
	
    
getGeneAgeEssentialityTable = function(geneAgeList, essentialGenes) {
    
    geneAgeEssList = cbind(geneAgeList, GeneIsEssential = FALSE)
    geneAgeEssList = as.data.frame(geneAgeEssList, stringsAsFactors = FALSE)
    	for (i in 1:nrow(essentialGenes)) {
		if (i %% 1000 == 0) print(paste0("Checked ", i, "/", nrow(essentialGenes), " rows."))
		
		thisGene = essentialGenes[i,1]
		geneAgeEssList[geneAgeEssList$ensembl_gene_id == thisGene, "GeneIsEssential"] = TRUE
		}

	rownames(geneAgeEssList) = NULL
	return(geneAgeEssList)
    
    }

tableReducer = function(myGeneTADAgeTable, ageNumbersFile, seed = 0) {

	set.seed(seed) # comment this if you don't want a deterministic random distribution
	
	ageNumbers = as.data.frame(read.table(ageNumbersFile, header = TRUE, sep = "\t"))
	ageNumbers$age = as.character(ageNumbers$age)
	ageNumbers$ngenes = as.numeric(ageNumbers$ngenes)
	
	resultingGeneTAGAgeTable = data.frame()
	for (i in 1:nrow(ageNumbers)) {
		
		thisAge = ageNumbers[i, 1]
		thisAgeN = ageNumbers[i, 2]
		
		print(paste0("Randomly picking ", thisAgeN, " genes of age ", thisAge))
		
		mySubset = subset(myGeneTADAgeTable, gene_age == thisAge)
		if (nrow(mySubset) > 0 && thisAgeN > 0) {
			if (nrow(mySubset) > thisAgeN) {
				mySubset$runif = runif(nrow(mySubset))
				mySubset = mySubset[order(mySubset$runif),]
				mySubset = mySubset[1:thisAgeN, 1:ncol(myGeneTADAgeTable)]
				}

			resultingGeneTAGAgeTable = rbind(resultingGeneTAGAgeTable, mySubset)
			
			}
		
		}
	rownames(resultingGeneTAGAgeTable) = NULL
	return(resultingGeneTAGAgeTable)
	
	}

################################################

# making chromosomes:

# to get centromere positions for human:
# go to http://genome.ucsc.edu/cgi-bin/hgTables
# choose assembly:"Feb. 2009 (GRCh37/hg19)" group:"All Tables" -> table:gap
# click on "get output"
# telomeres are coincidental with Caelinn's data, so I assume we are using the correct assembly
# filter centromeres are as gap "type"

# to get centromere positions for mouse:
# go to http://genome.ucsc.edu/cgi-bin/hgTables
# genome: Mouse
# choose assembly:"Dec. 2011 (GRCm38/mm10)" group:"All Tables" -> table:gap
# click on "get output"
# telomeres are coincidental with Caelinn's data, so I assume we are using the correct assembly
# filter centromeres are as gap "type"

centromerePositions_hg19 = as.data.frame(read.table(header = TRUE, sep = ",", text = "
chrom,centromere
chr1,121535434
chr2,92326171
chr3,90504854
chr4,49660117
chr5,46405641
chr6,58830166
chr7,58054331
chr8,43838887
chr9,47367679
chrX,58632012
chrY,10104553
chr10,39254935
chr11,51644205
chr12,34856694
chr13,16000000
chr14,16000000
chr15,17000000
chr16,35335801
chr17,22263006
chr18,15460898
chr19,24681782
chr20,26369569
chr21,11288129
chr22,13000000
"))

# for some reason, the ucsc website does not provide
# the end telomere of chr17, so I paste it here from Caelinn's list
endTelomerePositions_hg19 = as.data.frame(read.table(header = TRUE, sep = ",", text = "
chrom,end
chr1,249240621
chr2,243189373
chr3,198012430
chr4,191144276
chr5,180905260
chr6,171105067
chr7,159128663
chr8,146354022
chr9,141203431
chrX,155260560
chrY,59363566
chr10,135524747
chr11,134996516
chr12,133841895
chr13,115159878
chr14,107339540
chr15,102521392
chr16,90344753
chr17,81195210
chr18,78067248
chr19,59118983
chr20,63015520
chr21,48119895
chr22,51294566
"))


# no centromere for mouse chrY, so I give it the same value as for all others
centromerePositions_mm10 = as.data.frame(read.table(header = TRUE, sep = ",", text = "
chrom,centromere
chr1,110000
chr2,110000
chr3,110000
chr4,110000
chr5,110000
chr6,110000
chr7,110000
chr8,110000
chr9,110000
chrX,110000
chrY,110000
chr10,110000
chr11,110000
chr12,110000
chr13,110000
chr14,110000
chr15,110000
chr16,110000
chr17,110000
chr18,110000
chr19,110000
"))

endTelomerePositions_mm10 = as.data.frame(read.table(header = TRUE, sep = ",", text = "
chrom,end
chr1,195371971
chr2,182013224
chr3,159939680
chr4,156408116
chr5,151734684
chr6,149636546
chr7,145341459
chr8,129301213
chr9,124495110
chrX,170931299
chrY,91644698
chr10,130594993
chr11,121982543
chr12,120029022
chr13,120321639
chr14,124802244
chr15,103943685
chr16,98107768
chr17,94887271
chr18,90602639
chr19,61331566
"))

chromosomes_hg19 = vector()
for (i in 1:nrow(centromerePositions_hg19)) {
	myChrom = paste0(centromerePositions_hg19[i,]$chrom)
	myCentromere = centromerePositions_hg19[i,]$centromere
	myEnd = subset(endTelomerePositions_hg19, chrom == myChrom)$end
	chromosomes_hg19 = rbind(chromosomes_hg19, c(chrom = myChrom, centromere = myCentromere, end = myEnd))
	}
	
chromosomes_mm10 = vector()
for (i in 1:nrow(centromerePositions_mm10)) {
	myChrom = paste0(centromerePositions_mm10[i,]$chrom)
	myCentromere = centromerePositions_mm10[i,]$centromere
	myEnd = subset(endTelomerePositions_mm10, chrom == myChrom)$end
	chromosomes_mm10 = rbind(chromosomes_mm10, c(chrom = myChrom, centromere = myCentromere, end = myEnd))
	}

################################################
################################################
################################################
################################################

replacementTableMouse = as.data.frame(read.table(header = TRUE, sep = ",", text = "
from,to
Glires,Eutheria
Sciurognathi,Rodentia
Murinae,Rodentia
Mus.musculus,Rodentia
"))

replacementTableMouse_GliresRodentia = as.data.frame(read.table(header = TRUE, sep = ",", text = "
from,to
Glires,GliresRodentia
Rodentia,GliresRodentia
Sciurognathi,GliresRodentia
Murinae,GliresRodentia
Mus.musculus,GliresRodentia
"))

replacementTableHuman = as.data.frame(read.table(header = TRUE, sep = ",", text = "
from,to
Simiiformes,Primates
Catarrhini,Primates
Hominoidea,Primates
Hominidae,Primates
HomoPanGorilla,Primates
HomoSapiens,Primates
"))

ageListHuman = as.data.frame(read.table(header = TRUE, text = "
age
FungiMetazoa
Bilateria
Chordata
Euteleostomi
Sarcopterygii
Tetrapoda
Amniota
Mammalia
Theria
Eutheria
Primates
"))

ageListMouse = as.data.frame(read.table(header = TRUE, text = "
age
FungiMetazoa
Bilateria
Chordata
Euteleostomi
Sarcopterygii
Tetrapoda
Amniota
Mammalia
Theria
Eutheria
GliresRodentia
"))