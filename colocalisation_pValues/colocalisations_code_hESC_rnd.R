##############################
########## hESC_rnd ##########
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

geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_rndAge = randomiseAllAges2(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged, "gene_age")

####### final table, 16sep, hESC, rndAge
myFinalList_16sep_hESC_rndAge = getAllPairData(geneTADAgeTableAgeOnlyNoDupsHuman_PrimatesMerged_rndAge, ageListHuman, totRandomisations = 15)
MLMatrix_16sep_hESC_rndAge = fromFinalListToSgnTable(myFinalList_16sep_hESC_rndAge)
write.table(MLMatrix_16sep_hESC_rndAge, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_16sep_hESC_rndAge.tsv", sep = "\t")