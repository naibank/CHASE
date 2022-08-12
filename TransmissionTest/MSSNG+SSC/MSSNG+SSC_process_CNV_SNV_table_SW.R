################################################################################################################################################################################
# 
# MSSNG+SSC_process_CNV_SNV_table_SW.R
# purpose: combine MSSNG+SSC CNV_SNV_tables and add CH event frequency column 
# input: MSSNG_CNV_SNV_table.tsv, SSC_proband_CNV_SNV_table.tsv, SSC_unaffected_sibling_CNV_SNV_table.tsv
#         *database*.CNVs.freq.1perc.tsv, *database*.SVs.freq.1perc.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data
# output: MSSNG.SSC_CNV_SNV_table_eventfreq.tsv, SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv
#             -> used as input for MSSNG+SSC_get_FisherTest_snvfreq_SW.R
# 
# notes: modified from MSSNG_process_CNV_SNV_table.R & SSC_process_CNV_SNV_table.R (Faraz)
#
##############################################################################################################################################################################

library(data.table)


##############################################################################################################################################################################
### MSSNG CNV_SNV_table
setwd("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data")

MSSNG.CNV_SNV_table <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/MSSNG_CNV_SNV_table.tsv", data.table = F)
MSSNG.CNV_SNV_table$UID <- paste(MSSNG.CNV_SNV_table$`#Sample`, MSSNG.CNV_SNV_table$`#id`, MSSNG.CNV_SNV_table$sample, sep='.')
MSSNG.CNV_SNV_table <- dplyr::distinct(MSSNG.CNV_SNV_table, UID, .keep_all=T)
MSSNG.CNV_SNV_table$CNVID <- paste(MSSNG.CNV_SNV_table$sample, MSSNG.CNV_SNV_table$CHROM, MSSNG.CNV_SNV_table$START,
                             MSSNG.CNV_SNV_table$END)

cnv <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/MSSNG.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/MSSNG.SVs.freq.1perc.tsv", data.table = F)
cnv <- cnv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
sv <- sv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
cnvs <- rbind(cnv, sv)
colnames(cnvs)[5] <- 'CNV_freq'
cnvs$CNVID <- paste(cnvs$sample, cnvs$CHROM, cnvs$START,
                    cnvs$END)
## add new freq columns
MSSNG.CNV_SNV_table$CNV_freq <- NA
MSSNG.CNV_SNV_table$event_freq <- NA

# MSSNG.CNV_SNV_table <- merge(MSSNG.CNV_SNV_table, cnvs, by.x='sample', by.y='sample', all.x=F,all.y=F)
## add CNV freq column
for (row in 1:nrow(MSSNG.CNV_SNV_table)){
  message(row)
  MSSNG.CNV_SNV_table[row,]$CNV_freq <- cnvs$CNV_freq[cnvs$CNVID == MSSNG.CNV_SNV_table[row,]$CNVID][1]
}

## convert CNV percent to frequency
MSSNG.CNV_SNV_table$CNV_freq <- MSSNG.CNV_SNV_table$CNV_freq/100 


## add event_freq col
MSSNG.CNV_SNV_table$event_freq <- MSSNG.CNV_SNV_table$freq_max * MSSNG.CNV_SNV_table$CNV_freq ###Bank

# write.table(MSSNG.CNV_SNV_table, "MSSNG_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)

print("done MSSNG")

##############################################################################################################################################################################
### SSC Proband Table 

SSC.proband_table <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/SSC_proband_CNV_SNV_table.tsv", data.table = F)
SSC.proband_table$UID <- paste(SSC.proband_table$`#Sample`, SSC.proband_table$`#id`, SSC.proband_table$sample, sep='.')
SSC.proband_table <- dplyr::distinct(SSC.proband_table, UID, .keep_all=T)
SSC.proband_table$CNVID <- paste(SSC.proband_table$sample, SSC.proband_table$CHROM, SSC.proband_table$START,
                             SSC.proband_table$END)

cnv <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/SSC.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/SSC.SVs.freq.1perc.tsv", data.table = F)
cnv <- cnv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
sv <- sv[,c('sample', 'CHROM', 'START', 'END', 'freq_max')]
cnvs <- rbind(cnv, sv)
colnames(cnvs)[5] <- 'CNV_freq'
cnvs$CNVID <- paste(cnvs$sample, cnvs$CHROM, cnvs$START,
                    cnvs$END)
## add new freq columns
SSC.proband_table$CNV_freq <- NA
SSC.proband_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(SSC.proband_table)){
  message(row)
  SSC.proband_table[row,]$CNV_freq <- cnvs$CNV_freq[cnvs$CNVID == SSC.proband_table[row,]$CNVID][1]
}

## convert CNV percent to frequency
SSC.proband_table$CNV_freq <- SSC.proband_table$CNV_freq/100 


## add event_freq col
SSC.proband_table$event_freq <- SSC.proband_table$freq_max * SSC.proband_table$CNV_freq ###Bank

# write.table(SSC.proband_table, "SSC_proband_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)

print("done SSC proband")
print("start combine")


##############################################################################################################################################################################
# combine MSSNG and SSC proband CNV_SNV_tables
MSSNG.SSC.CNV_SNV_table <- rbind(MSSNG.CNV_SNV_table, SSC.proband_table)

write.table(MSSNG.SSC.CNV_SNV_table, "./MSSNG.SSC_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)
print("done combine")


##############################################################################################################################################################################
### SSC unaffected sibling Table 

SSC.unaffectedsib_table <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/SSC_unaffected_sibling_CNV_SNV_table.tsv", data.table = F)
SSC.unaffectedsib_table$UID <- paste(SSC.unaffectedsib_table$`#Sample`, SSC.unaffectedsib_table$`#id`, SSC.unaffectedsib_table$sample, sep='.')
SSC.unaffectedsib_table <- dplyr::distinct(SSC.unaffectedsib_table, UID, .keep_all=T)
SSC.unaffectedsib_table$CNVID <- paste(SSC.unaffectedsib_table$sample, SSC.unaffectedsib_table$CHROM, 
                                   SSC.unaffectedsib_table$START, SSC.unaffectedsib_table$END)
## add new freq columns
SSC.unaffectedsib_table$CNV_freq <- NA
SSC.unaffectedsib_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(SSC.unaffectedsib_table)){
  message(row)
  SSC.unaffectedsib_table[row,]$CNV_freq <- cnvs$CNV_freq[cnvs$CNVID == SSC.unaffectedsib_table[row,]$CNVID][1]
}

## convert CNV percent to frequency
SSC.unaffectedsib_table$CNV_freq <- SSC.unaffectedsib_table$CNV_freq/100 


## add event_freq col
SSC.unaffectedsib_table$event_freq <- SSC.unaffectedsib_table$freq_max * SSC.unaffectedsib_table$CNV_freq ###Bank

write.table(SSC.unaffectedsib_table, "./SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)

print("done SSC US")