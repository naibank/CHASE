################################################################################################################################################################################
# 
# MGRB_getSamplesWithCRVs_SW.R
#
# purpose: output list of MGRB samples with clinically relevant variants, including
#           Genetic Syndromes, Pathogenic SVs, and ASD Genes
# input: CNVs.MGRB.freq_1percent.HQR.tsv
#           /hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MGRB/CNVs/
#         MSSNG_SSC_pathogenic_SVs.xlsx, TADA MSSNG+SPARK+ASC ASD gene list - TADA MSSNG+SPARK+ASC ASD gene list.tsv
#           /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MGRB/data_processing/data/genomicSyndrome/
# output: samples.with.crCNV.MGRB.txt
#           /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MGRB/data_processing/data
# 
# notes:
#
##############################################################################################################################################################################

library(data.table)
library(gwascat)
library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BiocGenerics)
library(liftOver)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MGRB/data_processing/data/genomicSyndrome")


## Import hg38 to hg19 chain
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

# Steps:
# Remove cases with large CNVs > 3Mb, make sure to remove aneuploids
# Remove subjects with genomic syndromes, note the direction of CNV, deletion and duplication.
# Remove subjects with deletion impacting a gene on the genes in the list of clinically relevant CNVs.

cnvs <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MGRB/CNVs/CNVs.MGRB.freq_1percent.HQR.tsv", data.table = F)

#################################################################################################################################################
### Genetic Syndromes ####
dt <- cnvs [, c("sample", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")]
dt.del <- dt[which(dt$SVTYPE == "DEL"),]
dt.dup <- dt[which(dt$SVTYPE == "DUP"),]

big <- unique(dt$sample[dt$SIZE > 3000000])

gsd <- read.delim("./Genomic.Syndromes.hg38.20200606.txt", stringsAsFactors = F)
gsd.del <- gsd[which(gsd$cnvFlag == "Deletion"),]
gsd.dup <- gsd[which(gsd$cnvFlag == "Duplication"),]

dt.del.g <- GRanges(dt.del$chrAnn, IRanges(dt.del$STARTAnn, dt.del$ENDAnn), "*") # CNVs deletions
dt.dup.g <- GRanges(dt.dup$chrAnn, IRanges(dt.dup$STARTAnn, dt.dup$ENDAnn), "*") # CNVs duplications

gsd.del.hg38.g <- GRanges(gsd.del$chr, IRanges(gsd.del$start, gsd.del$end), "*") # genomic syndromes del
gsd.dup.hg38.g <- GRanges(gsd.dup$chr, IRanges(gsd.dup$start, gsd.dup$end), "*") # genomic syndromes dup
## Convert gsd GRanges (hg38) to hg19:
# export.bed(gsd.del.hg38.g,con='gsd.del.hg38.granges.bed')
# export.bed(gsd.dup.hg38.g,con='gsd.dup.hg38.granges.bed')

gsd.del.g <- import("hglft_genome_gsd.del.hg19.granges.bed")
gsd.del.g <- gsd.del.g[which(gsd.del.g$score == 1),]

gsd.dup.g <- import("hglft_genome_gsd.dup.hg19.granges.bed")
gsd.dup.g <- gsd.dup.g[which(gsd.dup.g$score == 1),]

# gsd.del.g <- unlist(liftOver(gsd.del.hg38.g, ch))
# gsd.dup.g <- unlist(liftOver(gsd.dup.hg38.g, ch))

## olap for deletions
olap.del <- data.frame(findOverlaps(dt.del.g, gsd.del.g))
olap.del$overlap <- width(pintersect(dt.del.g[olap.del$queryHits], gsd.del.g[olap.del$subjectHits])) # find #bp overlap for each pari
olap.del$overlap <- olap.del$overlap/width(gsd.del.g[olap.del$subjectHits]) # calculate % overlap
olap.del <- olap.del[olap.del$overlap > 0.9, ] # only include overlap > 90%
# olap.del$type <- ifelse(dt$SVTYPE[olap.del$queryHits] == "DEL", "Deletion", "Duplication")
# olap.del <- olap.del[olap.del$type == gsd$cnvFlag[olap.del$subjectHits], ]
gd.del <- unique(dt.del$sample[olap.del$queryHits])

## olap for duplications
olap.dup <- data.frame(findOverlaps(dt.dup.g, gsd.dup.g))
olap.dup$overlap <- width(pintersect(dt.dup.g[olap.dup$queryHits], gsd.dup.g[olap.dup$subjectHits])) # find #bp overlap for each pari
olap.dup$overlap <- olap.dup$overlap/width(gsd.dup.g[olap.dup$subjectHits]) # calculate % overlap
olap.dup <- olap.dup[olap.dup$overlap > 0.9, ] # only include overlap > 90%
# olap.dup$type <- ifelse(dt$SVTYPE[olap.dup$queryHits] == "DEL", "Deletion", "Duplication")
# olap.dup <- olap.dup[olap.dup$type == gsd.dup$cnvFlag[olap.dup$subjectHits], ]
gd.dup <- unique(dt.dup$sample[olap.dup$queryHits])

gd <- unique(append(gd.del, gd.dup))


#################################################################################################################################################
### Pathogenic SVs ####
crv <- readxl::read_excel("./MSSNG_SSC_pathogenic_SVs.xlsx")
crv <- crv$Annotation # gene and cnv type
crv <- crv[grep("del", crv)]
crv <- gsub(" del", "", crv)

inCRV <- function(genes, crv){
 return(sum(genes %in% crv) > 0) 
}

del <- dt[dt$SVTYPE == "DEL" & dt$exon_symbol != "" & !is.na(dt$exon_symbol), ]
del.genes <- sapply(sapply(del$exon_symbol, strsplit, "\\|"), inCRV, crv)
crv <- unique(del$sample[del.genes])


#################################################################################################################################################
### ASD Genes ####
asd <- fread("./TADA MSSNG+SPARK+ASC ASD gene list - TADA MSSNG+SPARK+ASC ASD gene list.tsv", data.table=F)

asd <- asd$gene 

inasd <- function(genes, asd){
  return(sum(genes %in% asd) > 0) 
}

del <- dt[dt$SVTYPE == "DEL" & dt$exon_symbol != "" & !is.na(dt$exon_symbol), ]
del.genes <- sapply(sapply(del$exon_symbol, strsplit, "\\|"), inasd, asd)
asd <- unique(del$sample[del.genes])


explained <- unique(c(big, crv, gd, asd))
writeLines(explained, "../samples.with.crCNV.MGRB.txt")
