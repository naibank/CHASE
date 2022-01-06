library(data.table)

CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)

## group LOF mutations
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                      CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF"

## filter table
CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$event_freq < 0.005 &
                                       CNV_SNV_table$Relation == "proband" &
                                       CNV_SNV_table$effect_priority == "LoF"),]

genes <- unique(CNV_SNV_table$gene_symbol)

write.table(genes, "SSC_eventfreq0.005_proband_genes.txt", quote = FALSE, row.names = FALSE, col.names = F)