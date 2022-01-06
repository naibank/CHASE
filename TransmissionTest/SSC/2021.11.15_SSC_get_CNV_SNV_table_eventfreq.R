library(data.table)
library(ggplot2)
#########################
##### Proband Table #####
#########################

proband_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_noCRV.tsv", data.table = F)
# proband_table <- fread("../data/SSC_proband_CNV_SNV_table_noCRV.tsv", data.table = F)

## remove CNV DGVpercFreq_subjects_coverageStudies_50percRecipOverlap > 10
# CNV_SNV_table <- CNV_SNV_table[which(!CNV_SNV_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap > 10),]

proband_table <- proband_table[which(proband_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap <= 10 
                                     | is.na(proband_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap)), ]

## add new freq columns
proband_table$CNV_freq <- NA
proband_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(proband_table)){
  proband_table[row,]$CNV_freq <- max(proband_table[row,25:31], na.rm = T) 
}

## convert CNV percent to frequency
proband_table$CNV_freq <- proband_table$CNV_freq/100


## add event_freq col
# for (row in 1:nrow(proband_table)){
  proband_table$event_freq <- proband_table$freq_max * proband_table$CNV_freq
# }

write.table(proband_table, "SSC_proband_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)



####################################
##### Unaffected Sibling Table #####
####################################

unaffectedsib_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table_noCRV.tsv", data.table = F)
# unaffectedsib_table <- fread("../data/SSC_unaffected_sibling_CNV_SNV_table_noCRV.tsv", data.table = F)

## remove CNV DGVpercFreq_subjects_coverageStudies_50percRecipOverlap > 10
# CNV_SNV_table <- CNV_SNV_table[which(!CNV_SNV_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap > 10),]

unaffectedsib_table <- unaffectedsib_table[which(unaffectedsib_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap <= 10 
                                     | is.na(unaffectedsib_table$DGVpercFreq_subjects_coverageStudies_50percRecipOverlap)), ]

## add new freq columns
unaffectedsib_table$CNV_freq <- NA
unaffectedsib_table$event_freq <- NA

## add CNV freq column
for (row in 1:nrow(unaffectedsib_table)){
  unaffectedsib_table[row,]$CNV_freq <- max(unaffectedsib_table[row,25:31], na.rm = T) 
}

## convert CNV percent to frequency
unaffectedsib_table$CNV_freq <- unaffectedsib_table$CNV_freq/100


## add event_freq col
# for (row in 1:nrow(unaffectedsib_table)){
  unaffectedsib_table$event_freq <- unaffectedsib_table$freq_max * unaffectedsib_table$CNV_freq
# }

write.table(unaffectedsib_table, "SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)


