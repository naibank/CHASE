library(data.table)
library(dplyr)


snvs <- fread("SSC_SNVs.tsv", data.table=F)
message('Done loading snvs')
cnv <- fread("../SSC.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("../SSC.SVs.freq.1perc.tsv", data.table = F)
cnv <- cnv[which(cnv$cds_symbol != ""),]
sv <- sv[which(sv$cds_symbol != "" & sv$length <= 10000),]
message('Done loading cnvs and svs')

cnv.genes <- unique(strsplit(paste(cnv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
sv.genes <- unique(strsplit(paste(sv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
snvs <- snvs[which(snvs$gene_symbol %in% c(cnv.genes, sv.genes)), ]

snvs_count <- nrow(snvs)
cnv_count <- nrow(cnv)
sv_count <- nrow(sv)

table <- data.frame(SNVCount=snvs_count, CNVCount=cnv_count, SVCount=sv_count)

write.table(table,"SSC_CNV_SNV_counts.tsv", sep="\t", row.names=F, quote=F, col.names=T)
