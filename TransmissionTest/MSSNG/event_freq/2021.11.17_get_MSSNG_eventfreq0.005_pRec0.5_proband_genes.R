library(data.table)

CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/CNV_SNV_table_eventfreq.tsv", data.table = F)


lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$event_freq < 0.005 &
                                       CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > 0.5] &
                                       CNV_SNV_table$Relation == "proband or affected sibling" &
                                       CNV_SNV_table$effect_priority == "nonsynonymous SNV"),]

genes <- unique(CNV_SNV_table$gene_symbol)

write.table(genes, "MSSNG_eventfreq0.005_pRec0.5_proband_genes.txt", quote = FALSE, row.names = FALSE, col.names = F)