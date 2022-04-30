library(data.table)
library(ggplot2)
#########################
##### Proband Table #####
#########################

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

proband_table <- fread("SSC_proband_CNV_SNV_table.tsv", data.table = F)
proband_table$UID <- paste(proband_table$`#Sample`, proband_table$`#id`, sep='.')
proband_table <- dplyr::distinct(proband_table, UID, .keep_all=T)
proband_table$CNVID <- paste(proband_table$sample, proband_table$CHROM, proband_table$START,
                             proband_table$END)

cnv <- fread("SSC.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("SSC.SVs.freq.1perc.tsv", data.table = F)
cnv <- cnv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
sv <- sv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
cnvs <- rbind(cnv, sv)
colnames(cnvs)[5] <- 'CNV_freq'
cnvs$CNVID <- paste(cnvs$sample, cnvs$CHROM, cnvs$START,
                    cnvs$END)
## add new freq columns
proband_table$CNV_freq <- NA
proband_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(proband_table)){
  message(row)
  proband_table[row,]$CNV_freq <- cnvs$CNV_freq[proband_table[row,]$CNVID == cnvs$CNVID]
}

## convert CNV percent to frequency
proband_table$CNV_freq <- proband_table$CNV_freq/100 


## add event_freq col
proband_table$event_freq <- proband_table$freq_max * proband_table$CNV_freq ###Bank

write.table(proband_table, "SSC_proband_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


####################################
##### Unaffected Sibling Table #####
####################################

unaffectedsib_table <- fread("SSC_unaffected_sibling_CNV_SNV_table.tsv", data.table = F)
unaffectedsib_table$UID <- paste(unaffectedsib_table$`#Sample`, unaffectedsib_table$`#id`, sep='.')
unaffectedsib_table <- dplyr::distinct(unaffectedsib_table, UID, .keep_all=T)
unaffectedsib_table$CNVID <- paste(unaffectedsib_table$sample, unaffectedsib_table$CHROM, 
                                   unaffectedsib_table$START, unaffectedsib_table$END)
## add new freq columns
unaffectedsib_table$CNV_freq <- NA
unaffectedsib_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(unaffectedsib_table)){
  message(row)
  unaffectedsib_table[row,]$CNV_freq <- cnvs$CNV_freq[unaffectedsib_table[row,]$CNVID == cnvs$CNVID]
}

## convert CNV percent to frequency
unaffectedsib_table$CNV_freq <- unaffectedsib_table$CNV_freq/100 


## add event_freq col
unaffectedsib_table$event_freq <- unaffectedsib_table$freq_max * unaffectedsib_table$CNV_freq ###Bank

write.table(unaffectedsib_table, "SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)



