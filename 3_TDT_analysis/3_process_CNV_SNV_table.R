library(data.table)

setwd("path_to_working_directory")

CNV_SNV_table <- fread("path_to_output_from_step_2", data.table = F)
CNV_SNV_table$UID <- paste(CNV_SNV_table$`#Sample`, CNV_SNV_table$`#id`, CNV_SNV_table$sample, sep='.')
CNV_SNV_table <- dplyr::distinct(CNV_SNV_table, UID, .keep_all=T)
CNV_SNV_table$CNVID <- paste(CNV_SNV_table$sample, CNV_SNV_table$CHROM, CNV_SNV_table$START,
                             CNV_SNV_table$END)

cnvs <- fread("path_to_all_CNVs_and_SVs", data.table = F)
colnames(cnvs)
colnames(cnvs)[11] <- 'CNV_freq'
cnvs$CNVID <- paste(cnvs$sample, cnvs$chrAnn, cnvs$STARTAnn, cnvs$ENDAnn)
# add new freq columns
CNV_SNV_table$CNV_freq <- NA
CNV_SNV_table$event_freq <- NA

# add CNV freq column
for (row in 1:nrow(CNV_SNV_table)){
  message(row)
  CNV_SNV_table[row,]$CNV_freq <- cnvs$CNV_freq[cnvs$CNVID == CNV_SNV_table[row,]$CNVID][1]
}

# convert CNV percent to frequency
CNV_SNV_table$CNV_freq <- CNV_SNV_table$CNV_freq/100 

# add event_freq col
CNV_SNV_table$event_freq <- CNV_SNV_table$freq_max * CNV_SNV_table$CNV_freq 

write.table(CNV_SNV_table, "CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


