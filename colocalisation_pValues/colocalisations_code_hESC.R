##############################
############ hESC ############
##############################
#
#

source("colocalisations_functions.R")

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