##############################
########## mESC_rnd ##########
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

geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_rndAge = randomiseAllAges2(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged, "gene_age")

####### final table, 05oct, mESC, rndAge
myFinalList_05oct_mESC_rndAge = getAllPairData(geneTADAgeTableAgeOnlyNoDupsMouse_RodentsMerged_rndAge, ageListMouse, totRandomisations = 15)
MLMatrix_05oct_mESC_rndAge = fromFinalListToSgnTable(myFinalList_05oct_mESC_rndAge)
write.table(MLMatrix_05oct_mESC_rndAge, "/mnt/sexreg/temporal/geneage/minusLogPValueSgn_table_05oct_mESC_rndAge.tsv", sep = "\t")