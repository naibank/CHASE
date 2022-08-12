library(data.table)
library(GenomicRanges)
# Steps:
# Remove cases with large CNVs > 3Mb, make sure to remove aneuploids
# Remove subjects with genomic syndromes, note the direction of CNV, deletion and duplication.
# Remove subjects with deletion impacting a gene on the genes in the list of clinically relevant CNVs.

cg <- fread("CNVs.CG/CNVs.MSSNG_CG.freq_10percent.HQR.IntFreq.tsv", data.table = F)
ilmn <- fread("CNVs.ILMN/CNVs.MSSNG_ILMN.freq_10percent.HQR.IntFreq.tsv", data.table = F)
ssc <- fread("SSC/CNVs.SSC.freq_10percent.HQR.IntFreq.tsv", data.table = F)

dt <- rbind(cg[, c("sample", "Family.ID", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")],
            ilmn[, c("sample", "Family.ID", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")],
            ssc[, c("sample", "Family.ID", "chrAnn", "STARTAnn", "ENDAnn", "SIZE", "SVTYPE", "exon_symbol", "exon_egID")])

big <- unique(dt$sample[dt$SIZE > 3000000])

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

explained <- unique(c(big, crv, gd))
writeLines(explained, "samples.with.crCNV.txt")
