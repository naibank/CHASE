################################################################################################################################################################################
# 
# SPARK_getSamplesWithCRVs_SW.R
#
# purpose: output list of SPARK samples with clinically relevant variants, including
#           Genetic Syndromes, Pathogenic SVs, and ASD Genes
# input: CNVs.SPARK_WGS_*1-3*.freq_1percent.HQR.tsv
#           /hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_1/CNVs/
#         Genomic.Syndromes.hg38.20200606.txt, MSSNG_SSC_pathogenic_SVs.xlsx, TADA MSSNG+SPARK+ASC ASD gene list - TADA MSSNG+SPARK+ASC ASD gene list.tsv
#           /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/SPARK/data_processing/data/genomicSyndrome/
# output: samples.with.crCNV.SPARK.txt
#           /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/SPARK/data_processing/data
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

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/SPARK/data_processing/data")

# Steps:
# Remove cases with large CNVs > 3Mb, make sure to remove aneuploids
# Remove subjects with genomic syndromes, note the direction of CNV, deletion and duplication.
# Remove subjects with deletion impacting a gene on the genes in the list of clinically relevant CNVs.

wgs1 <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_1/CNVs/CNVs.SPARK_WGS_1.freq_1percent.HQR.tsv", data.table = F)
wgs2 <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_2/CNVs/CNVs.SPARK_WGS_2.freq_1percent.HQR.tsv", data.table = F)
wgs3 <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_3/CNVs/CNVs.SPARK_WGS_3.freq_1percent.HQR.tsv", data.table = F)

dt <- rbind(wgs1[, c("sample", "Family ID", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")],
            wgs2[, c("sample", "Family ID", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")],
            wgs3[, c("sample", "Family ID", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")])

big <- unique(dt$sample[dt$SIZE > 3000000])

#################################################################################################################################################
### Genetic Syndromes ####
gsd <- read.delim("genomicSyndrome/Genomic.Syndromes.hg38.20200606.txt", stringsAsFactors = F)

dt.g <- GRanges(dt$chrAnn, IRanges(dt$STARTAnn, dt$ENDAnn), "*")
gsd.g <- GRanges(gsd$chr, IRanges(gsd$start, gsd$end), "*")
olap <- data.frame(findOverlaps(dt.g, gsd.g))
olap$overlap <- width(pintersect(dt.g[olap$queryHits], gsd.g[olap$subjectHits]))
olap$overlap <- olap$overlap/width(gsd.g[olap$subjectHits])
olap <- olap[olap$overlap > 0.9, ]
olap$type <- ifelse(dt$SVTYPE[olap$queryHits] == "DEL", "Deletion", "Duplication")
olap <- olap[olap$type == gsd$cnvFlag[olap$subjectHits], ]
gd <- unique(dt$sample[olap$queryHits])


#################################################################################################################################################
### Pathogenic SVs ####
crv <- readxl::read_excel("genomicSyndrome/MSSNG_SSC_pathogenic_SVs.xlsx")
crv <- crv$Annotation
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
asd <- fread("genomicSyndrome/TADA MSSNG+SPARK+ASC ASD gene list - TADA MSSNG+SPARK+ASC ASD gene list.tsv", data.table=F)

asd <- asd$gene 

inasd <- function(genes, asd){
  return(sum(genes %in% asd) > 0) 
}

del <- dt[dt$SVTYPE == "DEL" & dt$exon_symbol != "" & !is.na(dt$exon_symbol), ]
del.genes <- sapply(sapply(del$exon_symbol, strsplit, "\\|"), inasd, asd)
asd <- unique(del$sample[del.genes])


explained <- unique(c(big, crv, gd, asd))
 writeLines(explained, "./samples.with.crCNV.SPARK.txt")
