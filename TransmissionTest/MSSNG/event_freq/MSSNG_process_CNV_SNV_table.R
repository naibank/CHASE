library(data.table)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

CNV_SNV_table <- fread("MSSNG_CNV_SNV_table.tsv", data.table = F)
CNV_SNV_table$UID <- paste(CNV_SNV_table$`#Sample`, CNV_SNV_table$`#id`, sep='.')
CNV_SNV_table <- dplyr::distinct(CNV_SNV_table, UID, .keep_all=T)
CNV_SNV_table$CNVID <- paste(CNV_SNV_table$sample, CNV_SNV_table$CHROM, CNV_SNV_table$START,
                             CNV_SNV_table$END)

cnv <- fread("MSSNG.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("MSSNG.SVs.freq.1perc.tsv", data.table = F)
cnv <- cnv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
sv <- sv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
cnvs <- rbind(cnv, sv)
colnames(cnvs)[5] <- 'CNV_freq'
cnvs$CNVID <- paste(cnvs$sample, cnvs$CHROM, cnvs$START,
                    cnvs$END)
## add new freq columns
CNV_SNV_table$CNV_freq <- NA
CNV_SNV_table$event_freq <- NA

# CNV_SNV_table <- merge(CNV_SNV_table, cnvs, by.x='sample', by.y='sample', all.x=F,all.y=F)
## add CNV freq column
for (row in 1:nrow(CNV_SNV_table)){
  message(row)
  CNV_SNV_table[row,]$CNV_freq <- cnvs$CNV_freq[CNV_SNV_table[row,]$CNVID == cnvs$CNVID]
}

## convert CNV percent to frequency
CNV_SNV_table$CNV_freq <- CNV_SNV_table$CNV_freq/100 


## add event_freq col
CNV_SNV_table$event_freq <- CNV_SNV_table$freq_max * CNV_SNV_table$CNV_freq ###Bank

write.table(CNV_SNV_table, "MSSNG_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


