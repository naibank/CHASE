library(data.table)
load("gsMain_noLoFIn.Rdata")
meta <- read.delim("SSC_metadata.tsv", stringsAsFactors = F)
meta$Relation <- gsub("other sibling", "unaffected sibling", meta$Relation)
meta$Relation <- gsub("father|mother", "parent", meta$Relation)

proband_table <- fread("SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
unaffsib_table <- fread("SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)


lof <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

table.out <- data.frame()
  
for (gset in c("FMR1_Targets_Ascano", "PhMm_Aggr_CardvascMuscle_all")) {
  
  ## get target probands 
  target_proband <- proband_table[which(proband_table$event_freq < 0.0001 &
                                      proband_table$gene_symbol %in% lof$gene[lof$pRec > 0.9] &
                                      proband_table$effect_priority == "nonsynonymous SNV" &
                                      proband_table$entrez_id %in% gsMain[[gset]] & 
                                        proband_table$`#Sample` %in% meta$Sample.ID[meta$Relation == "proband"]),]
  
  ## get list of family IDs of probands in gene set
  target_famID <- meta$Family.ID[which(meta$Sample.ID %in% target_proband$`#Sample`)]
  
  ## get list of unaffected sibling IDs related to target proband
  related_unaffdsib <- meta$Sample.ID[which(meta$Family.ID %in% target_famID &
                                             meta$Relation == "unaffected sibling")]
  
  ## get target unaffected sibs (significant in target gene set)
  target_unaffsib <- unaffsib_table[which(unaffsib_table$event_freq < 0.0001 &
                                        unaffsib_table$gene_symbol %in% lof$gene[lof$pRec > 0.9] &
                                        unaffsib_table$effect_priority == "nonsynonymous SNV" &
                                        unaffsib_table$entrez_id %in% gsMain[[gset]] & 
                                          unaffsib_table$`#Sample` %in% meta$Sample.ID[meta$Relation == "unaffected sibling"]),]
  
  ## get target unaffected sibs that are related to target proband
  related_target_unaffsib <- sum(target_unaffsib$`#Sample` %in% related_unaffdsib)
  
  ## make table
  table.out <- rbind(table.out, 
                     data.frame("gene_set" = gset,
                                "target_proband" = nrow(target_proband),
                                "target_unaffected_sibling" = nrow(target_unaffsib),
                                "related_target_unaffected_slibing" = length(related_target_unaffsib)))

}

write.table(table.out, "SSC_unaffectedsib_relation", sep="\t", row.names=F, quote=F, col.names=T)

