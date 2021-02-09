#############################
######### mESC files ########
#############################

source("colocalisations_functions.R")

# I have checked that mm10 is equivalent to the feb2014 build also in mouse:
# see http://www.ensembl.org/Help/ArchiveRedirect?src=http%3A%2F%2Fjul2012.archive.ensembl.org%2F
library("biomaRt")
ensemblMouse = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='mmusculus_gene_ensembl')
allBioMartMouse = getBM(mart = ensemblMouse, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

thisFolder = "/mnt/sexreg/temporal/geneage/"

mESCTADs = read.table("/mnt/sexreg/ChromHMM/TADs/tsv/mESCTADs.tsv")
mESCTADs$chr = substr(mESCTADs$chr, 4, length(mESCTADs$chr) - 3)

sampleAgeFileMouse = paste0(thisFolder, "MiceAges_MTH.tsv")
# from Caelinn's email of 16mar-2020, file "MiceAges.txt", changing only Fungi-Metazoa_group to FungiMetazoa

geneAgeListMousePre = as.data.frame(read.table(sampleAgeFileMouse, header = TRUE, sep = "\t"))
geneAgeListMouse = as.data.frame(cbind(ensembl_gene_id = as.character(geneAgeListMousePre$GeneID), HUGO_symbol = NA, seqnames = NA, start = NA, end = NA, GeneAge = as.character(geneAgeListMousePre$GeneAge)))


newTADListMouse = getTADGenesTable(allBioMartMouse, mESCTADs)
geneTADAgeTableMouse = mergeGeneTADAge(newTADListMouse, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse = subset(geneTADAgeTableMouse, gene_age != "No Age Provided")
# ageList = getagelist(geneTADAgeTableAgeOnlyMouse) # but this has no order, so better use....

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

geneTADAgeTableAgeOnlyNoDupsMouse = geneTADAgeTableAgeOnlyMouse[!duplicated(geneTADAgeTableAgeOnlyMouse[, "ensembl_gene_id"]),]

geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse, 11, replacementTableMouse)

####### final table, 16sep, mESC, real
myFinalList_16sep_mESC_real = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged, ageListMouse, totRandomisations = 15)
MLMatrix_16sep_mESC_real = fromFinalListToSgnTable(myFinalList_16sep_mESC_real)
write.table(MLMatrix_16sep_mESC_real, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_mESC_real.tsv", sep = "\t")

MLMatrix_16sep_mESC_NH = fromFinalListToSgnTable(myFinalList_16sep_mESC_real, getHugeTADdataInstead = TRUE)
write.table(MLMatrix_16sep_mESC_NH, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_mESC_NH.tsv", sep = "\t")

ageNumbers_real = ngeneInAge(ageListMouse, geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged)
write.table(ageNumbers_real, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_mESC_real.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

ageNumbers_NH = ngeneInAge(ageListMouse, makeAllOneHugeTAD(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged))
write.table(ageNumbers_NH, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_mESC_NH.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


####### final table, 16sep, mESC, iTR = inter TAD removed
mESCTADs_iTR = interTADRemover(mESCTADs, chromosomes_mm10)
newTADListMouse_iTR = getTADGenesTable(allBioMartMouse, mESCTADs_iTR)
geneTADAgeTableMouse_iTR = mergeGeneTADAge(newTADListMouse_iTR, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse_iTR = subset(geneTADAgeTableMouse_iTR, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDupsMouse_iTR = geneTADAgeTableAgeOnlyMouse_iTR[!duplicated(geneTADAgeTableAgeOnlyMouse_iTR[, "ensembl_gene_id"]),]
geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_iTR = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse_iTR, 11, replacementTableMouse)

myFinalList_16sep_mESC_iTR = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_iTR, ageListMouse, totRandomisations = 15)
MLMatrix_16sep_mESC_iTR = fromFinalListToSgnTable(myFinalList_16sep_mESC_iTR)
write.table(MLMatrix_16sep_mESC_iTR, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_mESC_iTR.tsv", sep = "\t")

ageNumbers_16sep_mESC_iTR = ngeneInAge(ageListMouse, geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_iTR)
write.table(ageNumbers_16sep_mESC_iTR, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_mESC_iTR.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####### final table, 16sep, mESC, fused = fused TADs
mESCTADs_fused = fusedTAD(mESCTADs, chromosomes_mm10)
newTADListMouse_fused = getTADGenesTable(allBioMartMouse, mESCTADs_fused)
geneTADAgeTableMouse_fused = mergeGeneTADAge(newTADListMouse_fused, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse_fused = subset(geneTADAgeTableMouse_fused, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDupsMouse_fused = geneTADAgeTableAgeOnlyMouse_fused[!duplicated(geneTADAgeTableAgeOnlyMouse_fused[, "ensembl_gene_id"]),]
geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_fused = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse_fused, 11, replacementTableMouse)

myFinalList_16sep_mESC_fused = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_fused, ageListMouse, totRandomisations = 15)
MLMatrix_16sep_mESC_fused = fromFinalListToSgnTable(myFinalList_16sep_mESC_fused)
write.table(MLMatrix_16sep_mESC_fused, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_mESC_fused.tsv", sep = "\t")

ageNumbers_16sep_mESC_fused = ngeneInAge(ageListMouse, geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_fused)
write.table(ageNumbers_16sep_mESC_fused, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_mESC_fused.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#######################
#######################
####################### MOUSE ZONE GliresRodentia
#######################
#######################

# I have checked that mm10 is equivalent to the feb2014 build also in mouse:
# see http://www.ensembl.org/Help/ArchiveRedirect?src=http%3A%2F%2Fjul2012.archive.ensembl.org%2F
library("biomaRt")
ensemblMouse = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='mmusculus_gene_ensembl')
allBioMartMouse = getBM(mart = ensemblMouse, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

thisFolder = "/mnt/sexreg/temporal/geneage/"

mESCTADs = read.table("/mnt/sexreg/ChromHMM/TADs/tsv/mESCTADs.tsv")
mESCTADs$chr = substr(mESCTADs$chr, 4, length(mESCTADs$chr) - 3)

sampleAgeFileMouse = paste0(thisFolder, "MiceAges_MTH.tsv")
# from Caelinn's email of 16mar-2020, file "MiceAges.txt", changing only Fungi-Metazoa_group to FungiMetazoa

geneAgeListMousePre = as.data.frame(read.table(sampleAgeFileMouse, header = TRUE, sep = "\t"))
geneAgeListMouse = as.data.frame(cbind(ensembl_gene_id = as.character(geneAgeListMousePre$GeneID), HUGO_symbol = NA, seqnames = NA, start = NA, end = NA, GeneAge = as.character(geneAgeListMousePre$GeneAge)))


newTADListMouse = getTADGenesTable(allBioMartMouse, mESCTADs)
geneTADAgeTableMouse = mergeGeneTADAge(newTADListMouse, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse = subset(geneTADAgeTableMouse, gene_age != "No Age Provided")
# ageList = getagelist(geneTADAgeTableAgeOnlyMouse) # but this has no order, so better use....

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

geneTADAgeTableAgeOnlyNoDupsMouse = geneTADAgeTableAgeOnlyMouse[!duplicated(geneTADAgeTableAgeOnlyMouse[, "ensembl_gene_id"]),]

geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse, 11, replacementTableMouse_GliresRodentia)

####### final table, 05oct, mESC, real
myFinalList_05oct_mESC_real = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged, ageListMouse, totRandomisations = 15)
MLMatrix_05oct_mESC_real = fromFinalListToSgnTable(myFinalList_05oct_mESC_real)
write.table(MLMatrix_05oct_mESC_real, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_05oct_mESC_real.tsv", sep = "\t")

MLMatrix_05oct_mESC_NH = fromFinalListToSgnTable(myFinalList_05oct_mESC_real, getHugeTADdataInstead = TRUE)
write.table(MLMatrix_05oct_mESC_NH, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_05oct_mESC_NH.tsv", sep = "\t")

ageNumbers_real = ngeneInAge(ageListMouse, geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged)
write.table(ageNumbers_real, "/mnt/sexreg/temporal/geneage/ageNumbers_05oct_mESC_real.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

ageNumbers_NH = ngeneInAge(ageListMouse, makeAllOneHugeTAD(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged))
write.table(ageNumbers_NH, "/mnt/sexreg/temporal/geneage/ageNumbers_05oct_mESC_NH.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


####### final table, 05oct, mESC, iTR = inter TAD removed
mESCTADs_iTR = interTADRemover(mESCTADs, chromosomes_mm10)
newTADListMouse_iTR = getTADGenesTable(allBioMartMouse, mESCTADs_iTR)
geneTADAgeTableMouse_iTR = mergeGeneTADAge(newTADListMouse_iTR, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse_iTR = subset(geneTADAgeTableMouse_iTR, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDupsMouse_iTR = geneTADAgeTableAgeOnlyMouse_iTR[!duplicated(geneTADAgeTableAgeOnlyMouse_iTR[, "ensembl_gene_id"]),]
geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_iTR = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse_iTR, 11, replacementTableMouse_GliresRodentia)

myFinalList_05oct_mESC_iTR = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_iTR, ageListMouse, totRandomisations = 15)
MLMatrix_05oct_mESC_iTR = fromFinalListToSgnTable(myFinalList_05oct_mESC_iTR)
write.table(MLMatrix_05oct_mESC_iTR, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_05oct_mESC_iTR.tsv", sep = "\t")

ageNumbers_05oct_mESC_iTR = ngeneInAge(ageListMouse, geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_iTR)
write.table(ageNumbers_05oct_mESC_iTR, "/mnt/sexreg/temporal/geneage/ageNumbers_05oct_mESC_iTR.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####### final table, 05oct, mESC, fused = fused TADs
mESCTADs_fused = fusedTAD(mESCTADs, chromosomes_mm10)
newTADListMouse_fused = getTADGenesTable(allBioMartMouse, mESCTADs_fused)
geneTADAgeTableMouse_fused = mergeGeneTADAge(newTADListMouse_fused, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse_fused = subset(geneTADAgeTableMouse_fused, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDupsMouse_fused = geneTADAgeTableAgeOnlyMouse_fused[!duplicated(geneTADAgeTableAgeOnlyMouse_fused[, "ensembl_gene_id"]),]
geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_fused = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse_fused, 11, replacementTableMouse_GliresRodentia)

myFinalList_05oct_mESC_fused = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_fused, ageListMouse, totRandomisations = 15)
MLMatrix_05oct_mESC_fused = fromFinalListToSgnTable(myFinalList_05oct_mESC_fused)
write.table(MLMatrix_05oct_mESC_fused, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_05oct_mESC_fused.tsv", sep = "\t")

ageNumbers_05oct_mESC_fused = ngeneInAge(ageListMouse, geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_fused)
write.table(ageNumbers_05oct_mESC_fused, "/mnt/sexreg/temporal/geneage/ageNumbers_05oct_mESC_fused.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#######################
#######################
#######################
####################### hESC ZONE
#######################
#######################

hESCTADs = read.table("/mnt/sexreg/ChromHMM/TADs/tsv/hESCTADs.tsv")
hESCTADs$chr = substr(hESCTADs$chr, 4, length(hESCTADs$chr) - 3)


library(ggplot2)
library("biomaRt")
ensembl75 = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='hsapiens_gene_ensembl')
allBioMartHuman = getBM(mart = ensembl75, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

thisFolder = "/mnt/sexreg/temporal/geneage/"
sampleAgeFile = paste0(thisFolder, "AgesWithChrX.tsv") # converted from Caelinn's file in /data
geneAgeListHuman = as.data.frame(read.table(sampleAgeFile, header = TRUE, sep = "\t"))
colnames(geneAgeListHuman) = c("ensembl_gene_id", colnames(geneAgeListHuman)[2:6])
geneAgeListHuman = geneAgeListHuman[,1:6]

newTADList = getTADGenesTable(allBioMartHuman, hESCTADs)
geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeListHuman)
geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]

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

geneTADAgeTableAgeOnlyHuman = geneTADAgeTableAgeOnly

geneTADAgeTableAgeOnlyNoDupsHuman = geneTADAgeTableAgeOnlyHuman[!duplicated(geneTADAgeTableAgeOnlyHuman[, "ensembl_gene_id"]),]

geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged = replaceAges(geneTADAgeTableAgeOnlyNoDupsHuman, 11, replacementTableHuman)

####### final table, 16sep, hESC, real
myFinalList_16sep_hESC_real = getAllPairData(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged, ageListHuman, totRandomisations = 15)
MLMatrix_16sep_hESC_real = fromFinalListToSgnTable(myFinalList_16sep_hESC_real)
write.table(MLMatrix_16sep_hESC_real, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_hESC_real.tsv", sep = "\t")

MLMatrix_16sep_hESC_NH = fromFinalListToSgnTable(myFinalList_16sep_hESC_real, getHugeTADdataInstead = TRUE)
write.table(MLMatrix_16sep_hESC_NH, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_hESC_NH.tsv", sep = "\t")

ageNumbers_real = ngeneInAge(ageListHuman, geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged)
write.table(ageNumbers_real, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_hESC_real.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

ageNumbers_NH = ngeneInAge(ageListHuman, makeAllOneHugeTAD(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged))
write.table(ageNumbers_NH, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_hESC_NH.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


####### final table, 16sep, hESC, iTR = inter TAD removed
hESCTADs_iTR = interTADRemover(hESCTADs, chromosomes_hg19)
newTADListHuman_iTR = getTADGenesTable(allBioMartHuman, hESCTADs_iTR)
geneTADAgeTableHuman_iTR = mergeGeneTADAge(newTADListHuman_iTR, geneAgeListHuman)
geneTADAgeTableAgeOnlyHuman_iTR = subset(geneTADAgeTableHuman_iTR, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDupsHuman_iTR = geneTADAgeTableAgeOnlyHuman_iTR[!duplicated(geneTADAgeTableAgeOnlyHuman_iTR[, "ensembl_gene_id"]),]
geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_iTR = replaceAges(geneTADAgeTableAgeOnlyNoDupsHuman_iTR, 11, replacementTableHuman)

myFinalList_16sep_hESC_iTR = getAllPairData(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_iTR, ageListHuman, totRandomisations = 15)
MLMatrix_16sep_hESC_iTR = fromFinalListToSgnTable(myFinalList_16sep_hESC_iTR)
write.table(MLMatrix_16sep_hESC_iTR, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_hESC_iTR.tsv", sep = "\t")

ageNumbers_16sep_hESC_iTR = ngeneInAge(ageListHuman, geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_iTR)
write.table(ageNumbers_16sep_hESC_iTR, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_hESC_iTR.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####### final table, 16sep, hESC, fused = fused TADs
hESCTADs_fused = fusedTAD(hESCTADs, chromosomes_hg19)
newTADListHuman_fused = getTADGenesTable(allBioMartHuman, hESCTADs_fused)
geneTADAgeTableHuman_fused = mergeGeneTADAge(newTADListHuman_fused, geneAgeListHuman)
geneTADAgeTableAgeOnlyHuman_fused = subset(geneTADAgeTableHuman_fused, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDupsHuman_fused = geneTADAgeTableAgeOnlyHuman_fused[!duplicated(geneTADAgeTableAgeOnlyHuman_fused[, "ensembl_gene_id"]),]
geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_fused = replaceAges(geneTADAgeTableAgeOnlyNoDupsHuman_fused, 11, replacementTableHuman)

myFinalList_16sep_hESC_fused = getAllPairData(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_fused, ageListHuman, totRandomisations = 15)
MLMatrix_16sep_hESC_fused = fromFinalListToSgnTable(myFinalList_16sep_hESC_fused)
write.table(MLMatrix_16sep_hESC_fused, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_hESC_fused.tsv", sep = "\t")

ageNumbers_16sep_hESC_fused = ngeneInAge(ageListHuman, geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_fused)
write.table(ageNumbers_16sep_hESC_fused, "/mnt/sexreg/temporal/geneage/ageNumbers_16sep_hESC_fused.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#######################
#######################
#######################
####################### hESC ZONE rndAge
#######################
#######################

hESCTADs = read.table("/mnt/sexreg/ChromHMM/TADs/tsv/hESCTADs.tsv")
hESCTADs$chr = substr(hESCTADs$chr, 4, length(hESCTADs$chr) - 3)


library(ggplot2)
library("biomaRt")
ensembl75 = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='hsapiens_gene_ensembl')
allBioMartHuman = getBM(mart = ensembl75, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

thisFolder = "/mnt/sexreg/temporal/geneage/"
sampleAgeFile = paste0(thisFolder, "AgesWithChrX.tsv") # converted from Caelinn's file in /data
geneAgeListHuman = as.data.frame(read.table(sampleAgeFile, header = TRUE, sep = "\t"))
colnames(geneAgeListHuman) = c("ensembl_gene_id", colnames(geneAgeListHuman)[2:6])
geneAgeListHuman = geneAgeListHuman[,1:6]

newTADList = getTADGenesTable(allBioMartHuman, hESCTADs)
geneTADAgeTable = mergeGeneTADAge(newTADList, geneAgeListHuman)
geneTADAgeTableAgeOnly = subset(geneTADAgeTable, gene_age != "No Age Provided")
geneTADAgeTableAgeOnlyNoDups = geneTADAgeTableAgeOnly[!duplicated(geneTADAgeTableAgeOnly[, "ensembl_gene_id"]),]

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

geneTADAgeTableAgeOnlyHuman = geneTADAgeTableAgeOnly

geneTADAgeTableAgeOnlyNoDupsHuman = geneTADAgeTableAgeOnlyHuman[!duplicated(geneTADAgeTableAgeOnlyHuman[, "ensembl_gene_id"]),]

geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged = replaceAges(geneTADAgeTableAgeOnlyNoDupsHuman, 11, replacementTableHuman)

geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_rndAge = randomiseAllAges2(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged, "gene_age")

####### final table, 16sep, hESC, rndAge
myFinalList_16sep_hESC_rndAge = getAllPairData(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_rndAge, ageListHuman, totRandomisations = 15)
MLMatrix_16sep_hESC_rndAge = fromFinalListToSgnTable(myFinalList_16sep_hESC_rndAge)
write.table(MLMatrix_16sep_hESC_rndAge, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_hESC_rndAge.tsv", sep = "\t")

#######################
#######################
####################### MOUSE ZONE GliresRodentia rndAge
#######################
#######################

# I have checked that mm10 is equivalent to the feb2014 build also in mouse:
# see http://www.ensembl.org/Help/ArchiveRedirect?src=http%3A%2F%2Fjul2012.archive.ensembl.org%2F
library("biomaRt")
ensemblMouse = useMart(host='feb2014.archive.ensembl.org',
						biomart='ENSEMBL_MART_ENSEMBL',
						dataset='mmusculus_gene_ensembl')
allBioMartMouse = getBM(mart = ensemblMouse, attributes = c("ensembl_gene_id", "wikigene_name", "chromosome_name", "start_position", "end_position", "gene_biotype"))

thisFolder = "/mnt/sexreg/temporal/geneage/"

mESCTADs = read.table("/mnt/sexreg/ChromHMM/TADs/tsv/mESCTADs.tsv")
mESCTADs$chr = substr(mESCTADs$chr, 4, length(mESCTADs$chr) - 3)

sampleAgeFileMouse = paste0(thisFolder, "MiceAges_MTH.tsv")
# from Caelinn's email of 16mar-2020, file "MiceAges.txt", changing only Fungi-Metazoa_group to FungiMetazoa

geneAgeListMousePre = as.data.frame(read.table(sampleAgeFileMouse, header = TRUE, sep = "\t"))
geneAgeListMouse = as.data.frame(cbind(ensembl_gene_id = as.character(geneAgeListMousePre$GeneID), HUGO_symbol = NA, seqnames = NA, start = NA, end = NA, GeneAge = as.character(geneAgeListMousePre$GeneAge)))


newTADListMouse = getTADGenesTable(allBioMartMouse, mESCTADs)
geneTADAgeTableMouse = mergeGeneTADAge(newTADListMouse, geneAgeListMouse)
geneTADAgeTableAgeOnlyMouse = subset(geneTADAgeTableMouse, gene_age != "No Age Provided")
# ageList = getagelist(geneTADAgeTableAgeOnlyMouse) # but this has no order, so better use....

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

geneTADAgeTableAgeOnlyNoDupsMouse = geneTADAgeTableAgeOnlyMouse[!duplicated(geneTADAgeTableAgeOnlyMouse[, "ensembl_gene_id"]),]

geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged = replaceAges(geneTADAgeTableAgeOnlyNoDupsMouse, 11, replacementTableMouse_GliresRodentia)

geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_rndAge = randomiseAllAges2(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged, "gene_age")

####### final table, 05oct, mESC, rndAge
myFinalList_05oct_mESC_rndAge = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_rndAge, ageListMouse, totRandomisations = 15)
MLMatrix_05oct_mESC_rndAge = fromFinalListToSgnTable(myFinalList_05oct_mESC_rndAge)
write.table(MLMatrix_05oct_mESC_rndAge, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_05oct_mESC_rndAge.tsv", sep = "\t")