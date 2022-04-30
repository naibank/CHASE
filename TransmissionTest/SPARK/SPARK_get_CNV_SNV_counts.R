library(data.table)
library(dplyr)


snvs <- rbind(fread("SPARKWGS1_SNVs.tsv", data.table=F),
              fread("SPARKWGS2_SNVs.tsv", data.table=F),
              fread("SPARKWGS3_SNVs.tsv", data.table=F))
message('Done loading snvs')
cnv <- rbind(fread("SPARKWGS1.CNVs.freq.1perc.tsv", data.table = F),
             fread("SPARKWGS2.CNVs.freq.1perc.tsv", data.table = F),
             fread("SPARKWGS3.CNVs.freq.1perc.tsv", data.table = F))
cnv <- cnv[which(cnv$cds_symbol != ""),]
message('Done loading cnvs')

cnv.genes <- unique(strsplit(paste(cnv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
snvs <- snvs[which(snvs$gene_symbol %in% c(cnv.genes)), ]

snvs_count <- nrow(snvs)
cnv_count <- nrow(cnv)

table <- data.frame(SNVCount=snvs_count, CNVCount=cnv_count, SVCount=NA)

write.table(table,"SPARK_CNV_SNV_counts.tsv", sep="\t", row.names=F, quote=F, col.names=T)
