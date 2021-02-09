##############################
############ mESC ############
##############################
#
# mESC merging GliresRodentia

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