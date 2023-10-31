ageEquivalencies_binary = as.data.frame(read.table(header = TRUE, sep = ",", text = "
age,maxMy,minMy,avgMy
FungiMetazoa,3,3,3
Bilateria,3,3,3
Chordata,3,3,3
Euteleostomi,3,3,3
Sarcopterygii,1,1,1
Tetrapoda,1,1,1
Amniota,1,1,1
Mammalia,1,1,1
Theria,1,1,1
Eutheria,1,1,1
Simiiformes,1,1,1
Catarrhini,1,1,1
Hominoidea,1,1,1
Hominidae,1,1,1
HomoPanGorilla,1,1,1
HomoSapiens,1,1,1
Glires,1,1,1
Rodentia,1,1,1
Sciurognathi,1,1,1
Murinae,1,1,1
Mus.musculus,1,1,1
"))

######################

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
		mySubset = subset(myTADs, myTADs$end > row$start_position & myTADs$start < row$end_position & substr(myTADs$chr, 4, 5) == row$chromosome_name)
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

getTADAgeList = function(myGeneTADAgeTable, ageEquival) {
	
	# using merge here is overkill, because I just want to bring the avgMy
	myGeneTADAgeTable$avgMy = ageEquival$avgMy[match(myGeneTADAgeTable$gene_age, ageEquival$age)]
	allTADs = as.data.frame(sort(as.character(unique(myGeneTADAgeTable$tad_id))))
	colnames(allTADs) = "tad_id"
	allTADs$tad_chr = NA
	allTADs$tad_start = NA
	allTADs$tad_end = NA
	allTADs$meanAge = NA
	allTADs$percentOld = NA

	for(i in 1:nrow(allTADs)) {
		allTADs[i, "meanAge"] = mean(subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"])[, "avgMy"])
		if (nrow(subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"])) > 0)
			allTADs[i, "percentOld"] =
				nrow(subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"] & avgMy == 3)) / nrow(subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"]))
		allTADs[i, "tad_chr"] = as.character(subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"])[1, "tad_chr"])
		allTADs[i, "tad_start"] = subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"])[1, "tad_start"]
		allTADs[i, "tad_end"] = subset(myGeneTADAgeTable, as.character(tad_id) == allTADs[i, "tad_id"])[1, "tad_end"]
		}

	allTADs = allTADs[with(allTADs, order(tad_chr, tad_start, tad_end)),]
	rownames(allTADs) = NULL
	
	return(allTADs)
	}

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

# function to calculate age ratios in stable or unstable TADs

get_ageratios = function(TADFile, ageFile, randomisations = 0, stableText = "stable") {

# sampleTADFile = "TADs_100kbBookendBoundaries_byStability_CTgt5_stable.bed"
# sampleAgeFile = "AgesWithChrX.tsv" # converted from Caelinn's file in /data

	TADList = as.data.frame(read.table(TADFile, header = TRUE))
	geneAgeList = as.data.frame(read.table(ageFile, header = TRUE, sep = "\t"))
	colnames(geneAgeList) = c("ensembl_gene_id", colnames(geneAgeList)[2:6])
	geneAgeList = geneAgeList[,1:6]

	newTADList = getTADGenesTable(allBioMartHuman, TADList)
	
	# randomisation == 0 means that it is not randomised, the other are randomised with seed = i
	for (i in 0:randomisations) {
		
		print(paste0("Working on loop ", i, "/", randomisations))
		
		filenameOut = paste0("TADs_percentOld_", stableText)
		mainText_base = paste0("Percent of old genes in TAD\nfor ", stableText, " TADs\n")
		
		geneAgeList_updated = geneAgeList
		filenameTable = paste0(filenameOut, ".txt")
		filenamePNG = paste0(filenameOut, ".png")
		mainText = mainText_base
		
		if (i > 0) {
			geneAgeList_updated = randomiseAllAges(geneAgeList, seed = i, ageColumn = "GeneAge", chrColumn = "seqnames")
			mainText = paste0(mainText, "randomisation with seed = ", i)
			rndText = sprintf(paste0("%0", floor(log10(randomisations)) + 1, "d"), i)
			filenameTable = paste0(filenameOut, "_rndseed=", rndText, ".txt")
			filenamePNG = paste0(filenameOut, "_rndseed=", rndText, ".png")
			}
		
		geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeList_updated)
		geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")
		geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]


		# get a list {chromosome, TAD id, average age of contained genes in My}
		TADAgeList = getTADAgeList(geneTADAgeTableAgeOnlyNoDups, ageEquivalencies_binary)

		write.table(TADAgeList, filenameTable, sep = "\t", quote = FALSE, row.names = FALSE)
		png(filenamePNG)
		hist(TADAgeList$percentOld, col = "salmon", main = mainText, xlim = c(-0.1, 1.1), breaks = (0:10)/10)
		dev.off()
		}
	
	}

######################


library("biomaRt")
ensembl75 = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='hsapiens_gene_ensembl')
allBioMartHuman = getBM(mart = ensembl75, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

setwd("./")

#
#
# calculating for stable TADs

sampleTADFile = "TADs_100kbBookendBoundaries_byStability_CTgt5_stable.bed"
sampleAgeFile = "AgesWithChrX.tsv" # converted from Caelinn's file in /data
TADList = as.data.frame(read.table(sampleTADFile, header = TRUE))
geneAgeList = as.data.frame(read.table(sampleAgeFile, header = TRUE, sep = "\t"))
colnames(geneAgeList) = c("ensembl_gene_id", colnames(geneAgeList)[2:6])
geneAgeList = geneAgeList[,1:6]

newTADList = getTADGenesTable(allBioMartHuman, TADList)
geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeList)
geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]


# get a list {chromosome, TAD id, average age of contained genes in My}
TADAgeList = getTADAgeList(geneTADAgeTableAgeOnlyNoDups, ageEquivalencies_binary)

write.table(TADAgeList, "TADs_meanAge_stable.txt", sep = "\t", quote = FALSE, row.names = FALSE)
png("TADs_meanAge_stable.png")
hist(TADAgeList$meanAge, col = "salmon", main = "Mean ages for stable TADs\n(young = 1, old = 3)")
dev.off()

write.table(TADAgeList, "TADs_percentOld_stable.txt", sep = "\t", quote = FALSE, row.names = FALSE)
png("TADs_percentOld_stable.png")
hist(TADAgeList$percentOld, col = "salmon", main = "Percent of old genes in TAD\nfor stable TADs\n", xlim = c(-0.1, 1.1), breaks = (0:10)/10)
dev.off()

#
#
# calculating for UNstable TADs

sampleTADFile = "TADs_100kbBookendBoundaries_byStability_CTgt5_unstable.bed"
sampleAgeFile = "AgesWithChrX.tsv" # converted from Caelinn's file in /data
TADList = as.data.frame(read.table(sampleTADFile, header = TRUE))
geneAgeList = as.data.frame(read.table(sampleAgeFile, header = TRUE, sep = "\t"))
colnames(geneAgeList) = c("ensembl_gene_id", colnames(geneAgeList)[2:6])
geneAgeList = geneAgeList[,1:6]

newTADList = getTADGenesTable(allBioMartHuman, TADList)
geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeList)
geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]


# get a list {chromosome, TAD id, average age of contained genes in My}
TADAgeList = getTADAgeList(geneTADAgeTableAgeOnlyNoDups, ageEquivalencies_binary)

write.table(TADAgeList, "TADs_meanAge_unstable.txt", sep = "\t", quote = FALSE, row.names = FALSE)
png("TADs_meanAge_unstable.png")
hist(TADAgeList$meanAge, col = "salmon", main = "Mean ages for unstable TADs\n(young = 1, old = 3)")
dev.off()

#
#
# calculating for UNstable TADs, second version (with TADs that are the same size in unstable as in stable)

sampleTADFile = "TADs_100kbBookendBoundaries_byStability_CTlt5_unstable_version2.bed"
sampleAgeFile = "AgesWithChrX.tsv" # converted from Caelinn's file in /data
TADList = as.data.frame(read.table(sampleTADFile, header = TRUE))
geneAgeList = as.data.frame(read.table(sampleAgeFile, header = TRUE, sep = "\t"))
colnames(geneAgeList) = c("ensembl_gene_id", colnames(geneAgeList)[2:6])
geneAgeList = geneAgeList[,1:6]

newTADList = getTADGenesTable(allBioMartHuman, TADList)
geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeList)
geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]


# get a list {chromosome, TAD id, average age of contained genes in My}
TADAgeList = getTADAgeList(geneTADAgeTableAgeOnlyNoDups, ageEquivalencies_binary)

write.table(TADAgeList, "TADs_meanAge_unstable_version2.txt", sep = "\t", quote = FALSE, row.names = FALSE)
png("TADs_meanAge_unstable_version2.png")
hist(TADAgeList$meanAge, col = "salmon", main = "Mean ages for unstable TADs\n[second definition]\n(young = 1, old = 3)")
dev.off()

write.table(TADAgeList, "TADs_percentOld_unstable_version2.txt", sep = "\t", quote = FALSE, row.names = FALSE)
png("TADs_percentOld_unstable_version2.png")
hist(TADAgeList$percentOld, col = "salmon", main = "Percent of old genes in TAD\nfor unstable TADs\n[second definition]", xlim = c(-0.1, 1.1), breaks = (0:10)/10)
dev.off()

######################

TADs_stable = as.data.frame(read.table("TADs_percentOld_stable.txt", header = TRUE, sep = "\t"))
TADs_unstable = as.data.frame(read.table("TADs_percentOld_unstable_version2.txt", header = TRUE, sep = "\t"))

######################

# USING FUNCTION

library("biomaRt")
ensembl75 = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='hsapiens_gene_ensembl')
allBioMartHuman = getBM(mart = ensembl75, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

setwd("./")

for (analysisType in c("CTgt5_stable", "CTlt5_unstable_version2")) {
	get_ageratios(	TADFile = paste0("TADs_100kbBookendBoundaries_byStability_", analysisType, ".bed"),
					ageFile = "AgesWithChrX.tsv",
					randomisations = 100,
					stableText = analysisType
					)
	}

######################

library(ggplot2)
setwd("./")
randomisations = 100

result = data.frame()
for (analysisType in c("CTgt5_stable", "CTlt5_unstable_version2")) {
	for (i in 0:randomisations) {
		filename = paste0("TADs_percentOld_", analysisType, ".txt")
		if (i > 0) filename = paste0("TADs_percentOld_", analysisType, "_rndseed=", sprintf(paste0("%0", floor(log10(randomisations)) + 1, "d"), i),".txt")
		myTable = read.table(filename, sep = "\t", header = TRUE)
		totTADs = nrow(myTable)
		totYoungTADs = nrow(subset(myTable, percentOld < 0.1))
		ratioYoungTADs = totYoungTADs / totTADs
		result = rbind(result, data.frame(file = filename, type = analysisType, randomised = i > 0, ratio_young = ratioYoungTADs))
		}
	}

p = ggplot(result) +
	geom_boxplot(aes(type, ratio_young), outlier.shape = NA) +
	geom_jitter(aes(type, ratio_young, colour = randomised),
		position = position_jitter(width = .15, height = -0.7), size = 2) +
	scale_y_continuous(limits = c(0, 0.20)) +
	scale_color_manual(values = c("darkred", "steelblue")) +
	labs(title = "Ratio of young TADs",
			subtitle = "TADs with < 10% of old genes,\ncomparing stable (left) and unstable TADs (right),\nobserved (red) vs 100 randomisations (blue)")
	
ggsave("boxplots_randomisations.png", plot = p, width = 8, height = 8, dpi = 100)
write.table(result, "boxplots_randomisations.txt", sep = "\t", quote = FALSE, row.names = FALSE)

######################

library(ggplot2)
setwd("./")
result = read.table("boxplots_randomisations.txt", sep = "\t", header = TRUE)

set.seed(0)
p = ggplot(result) +
	#geom_boxplot(aes(type, ratio_young), outlier.shape = NA) +
	geom_jitter(aes(type, ratio_young, colour = randomised),
		position = position_jitter(width = .15, height = -0.7), size = 2) +
	scale_y_continuous(limits = c(0, 0.20)) +
	scale_color_manual(values = c("darkred", "steelblue")) +
	labs(title = "Ratio of young TADs",
			subtitle = "TADs with < 10% of old genes,\ncomparing stable (left) and unstable TADs (right),\nobserved (red) vs 100 randomisations (blue)")

ggsave("boxplots_randomisations_jitter_only_0.png", plot = p, width = 8, height = 8, dpi = 100)