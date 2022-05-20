library(data.table)
library(ggplot2)
#########################
##### Proband Table #####
#########################

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/TT/")

Determine_CNV_MAF <- function(CNVs) {
  CNV_MAF <- pmax(CNVs$CGparentalPercFreq_50percRecipOverlap,
                  pmax(CNVs$cnvnIlmXParentalPercFreq_50percRecipOverlap,
                       pmax(CNVs$erdsIlmXParentalPercFreq_50percRecipOverlap,
                            pmax(CNVs$cnvnIlm2ParentalPercFreq_50percRecipOverlap,
                                 pmax(CNVs$erdsIlm2ParentalPercFreq_50percRecipOverlap,
                                      pmax(CNVs$otgCnvnPercFreq_50percRecipOverlap,
                                           pmax(CNVs$otgErdsPercFreq_50percRecipOverlap,
                                                pmax(CNVs$sscErdsPercFreq_50percRecipOverlap,
                                                     CNVs$sscCnvnPercFreq_50percRecipOverlap, na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T)/100
  return (CNV_MAF)
}

proband_table <- rbind(fread("./SPARKWGS1_proband_CNV_SNV_table.tsv", data.table = F),
                       fread("./SPARKWGS2_proband_CNV_SNV_table.tsv", data.table = F),
                       fread("./SPARKWGS3_proband_CNV_SNV_table.tsv", data.table = F))
proband_table$UID <- paste(proband_table$`#Sample`, proband_table$sample, proband_table$`#id`, sep='.')
# proband_table <- dplyr::distinct(proband_table, UID, .keep_all=T)
proband_table$CNVID <- paste(proband_table$sample, proband_table$CHROM, proband_table$START,
                             proband_table$END)

cnv <- rbind(fread("./SPARKWGS1.CNVs.freq.1perc.tsv", data.table = F),
             fread("./SPARKWGS2.CNVs.freq.1perc.tsv", data.table = F),
             fread("./SPARKWGS3.CNVs.freq.1perc.tsv", data.table = F))
cnv$freq_max <- Determine_CNV_MAF(cnv)
cnv <- cnv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
cnvs <- cnv
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


## add event_freq col
proband_table$event_freq <- proband_table$freq_max * proband_table$CNV_freq ###Bank

write.table(proband_table, "./SPARK_proband_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


####################################
##### Unaffected Sibling Table #####
####################################

unaffectedsib_table <- rbind(fread("./SPARKWGS1_unaffected_sibling_CNV_SNV_table.tsv", data.table = F),
                       fread("./SPARKWGS2_unaffected_sibling_CNV_SNV_table.tsv", data.table = F),
                       fread("./SPARKWGS3_unaffected_sibling_CNV_SNV_table.tsv", data.table = F))
unaffectedsib_table$UID <- paste(unaffectedsib_table$`#Sample`, unaffectedsib_table$sample, unaffectedsib_table$`#id`, sep='.')
# unaffectedsib_table <- dplyr::distinct(unaffectedsib_table, UID, .keep_all=T)
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

## add event_freq col
unaffectedsib_table$event_freq <- unaffectedsib_table$freq_max * unaffectedsib_table$CNV_freq ###Bank

write.table(unaffectedsib_table, "./SPARK_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


