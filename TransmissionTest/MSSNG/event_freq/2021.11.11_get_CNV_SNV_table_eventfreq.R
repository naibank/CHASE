library(data.table)

CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/2021.11.01 no CRVs/data/2021.11.01_CNV_SNV_table.tsv", data.table = F)
# CNV_SNV_table <- fread("../data/2021.11.01_CNV_SNV_table.tsv", data.table = F)[, -1]
# CNV_SNV_table$V1 <- NULL

## remove CNV DGVpercFreq_subjects_coverageStudies_50percRecipOverlap > 10
# CNV_SNV_table <- CNV_SNV_table[which(!CNV_SNV_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap > 10),]

CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap <= 10 
                                     | is.na(CNV_SNV_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap)), ]

## add new freq columns
CNV_SNV_table$CNV_freq <- NA
CNV_SNV_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(CNV_SNV_table)){
  CNV_SNV_table[row,]$CNV_freq <- max(CNV_SNV_table[row,25:31], na.rm = T) 
}

## convert CNV percent to frequency
CNV_SNV_table$CNV_freq <- CNV_SNV_table$CNV_freq/100 


## add event_freq col
CNV_SNV_table$event_freq <- CNV_SNV_table$freq_max * CNV_SNV_table$CNV_freq ###Bank

write.table(CNV_SNV_table, "CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


