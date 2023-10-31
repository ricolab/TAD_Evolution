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

chromosomes_hg19 = vector()
for (i in 1:nrow(centromerePositions_hg19)) {
	myChrom = paste0(centromerePositions_hg19[i,]$chrom)
	myCentromere = centromerePositions_hg19[i,]$centromere
	myEnd = subset(endTelomerePositions_hg19, chrom == myChrom)$end
	chromosomes_hg19 = rbind(chromosomes_hg19, c(chrom = myChrom, centromere = myCentromere, end = myEnd))
	}

chromosomes_hg19 = as.data.frame(chromosomes_hg19)

######################

setwd("./")

boundaries = as.data.frame(read.table("100kbBookendBoundaries_byStability.bed", header = TRUE))

getTADs_chr = function(mytable, chrlist) {
	outputtable = data.frame()
	chr = unique(mytable$chr)
	chrend = as.numeric(as.character(subset(chrlist, as.character(chrom) == as.character(chr))[,"end"]))
	for (i in 1:nrow(mytable)) {
		if (i == 1) start = 0
		else start = end
		end = round((mytable[i, "loc"] + mytable[i, "loc2"]) / 2)
		newrow = data.frame(chr = chr, start = start, end = end)
		outputtable = rbind(outputtable, newrow)
		}
	# add last TAD
	newrow = data.frame(chr = chr, start = end, end = chrend)
	outputtable = rbind(outputtable, newrow)
	return(outputtable)
	}

getTADs = function(boundaries_table, min_celltypes = 5, chrlist, greaterthan_or_equalto = TRUE) {
	if (greaterthan_or_equalto) boundaries_sub = subset(boundaries_table, counts >= min_celltypes)
	if (!greaterthan_or_equalto) boundaries_sub = subset(boundaries_table, counts < min_celltypes)
	result = data.frame()
	thesechr = as.character(unique(boundaries_sub$chr))
	for (thischr in thesechr) {
		subtable = subset(boundaries_sub, chr == thischr)
		subTADTable = getTADs_chr(subtable, chrlist)
		result = rbind(result, subTADTable)
		}
	return(result)
	}

finalTADs = getTADs(boundaries, min_celltypes = 5, chrlist = chromosomes_hg19)

write.table(finalTADs, file = "TADs_100kbBookendBoundaries_byStability_CTgt5.bed", sep = "\t", quote = FALSE, row.names = FALSE)

######################

# Creating the TADs that are unstable

finalTADs_unstable = getTADs(boundaries, min_celltypes = 5, chrlist = chromosomes_hg19, greaterthan_or_equalto = FALSE)
write.table(finalTADs_unstable, file = "TADs_100kbBookendBoundaries_byStability_CTlt5.bed", sep = "\t", quote = FALSE, row.names = FALSE)


######################

# Separating the different types of TADs

setwd("./")
boundaries = as.data.frame(read.table("100kbBookendBoundaries_byStability.bed", header = TRUE))
TADs = as.data.frame(read.table("TADs_100kbBookendBoundaries_byStability_CTgt5.bed", header = TRUE))

ctlimit = 5
TADs_stable = data.frame()
TADs_unstable = data.frame()
for (i in 1:nrow(TADs)) {
	TADStart = TADs[i, "start"]
	TADEnd = TADs[i, "end"]
	TADChr = TADs[i, "chr"]
	# check if there is a boundary with < ctlimit within the TAD
	TAD_is_stable = nrow(subset(boundaries, loc > TADStart & loc < TADEnd & loc2 > TADStart & loc2 < TADEnd & chr == TADChr & counts < ctlimit)) == 0
	if (TAD_is_stable) TADs_stable = rbind(TADs_stable, TADs[i,])
	else TADs_unstable = rbind(TADs_unstable, TADs[i,])
	}

write.table(TADs_stable, "TADs_100kbBookendBoundaries_byStability_CTgt5_stable.bed", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(TADs_unstable, "TADs_100kbBookendBoundaries_byStability_CTgt5_unstable.bed", sep = "\t", row.names = FALSE, quote = FALSE)

######################

# same for unstable

setwd("./")
boundaries = as.data.frame(read.table("100kbBookendBoundaries_byStability.bed", header = TRUE))
TADs = as.data.frame(read.table("TADs_100kbBookendBoundaries_byStability_CTlt5.bed", header = TRUE)) # note that is lt == «less than»

TADs_unstable = data.frame()
for (i in 1:nrow(TADs)) {
	TADStart = TADs[i, "start"]
	TADEnd = TADs[i, "end"]
	TADChr = TADs[i, "chr"]
	# check if there is a boundary within the TAD (regardless of ctlimit)
	TAD_is_empty = nrow(subset(boundaries, loc > TADStart & loc < TADEnd & loc2 > TADStart & loc2 < TADEnd & chr == TADChr)) == 0
	if (TAD_is_empty) TADs_unstable = rbind(TADs_unstable, TADs[i,])
	}

write.table(TADs_unstable, "TADs_100kbBookendBoundaries_byStability_CTlt5_unstable_version2.bed", sep = "\t", row.names = FALSE, quote = FALSE)