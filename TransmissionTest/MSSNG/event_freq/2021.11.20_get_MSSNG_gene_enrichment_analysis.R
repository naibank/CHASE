library(data.table)

CNV_SNV_table <- data.table::fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/CNV_SNV_table_eventfreq.tsv", data.table = F)

## group LOF variants in effect_priority
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                      CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank

test.out <- data.frame()

for (gsi in 1:length(gsMain)){
  for(freq in c(0.005, 0.001)){ # event freq cut-offs; 0.005 = no cut-off
    for(variant_type in c("LoF", "nonsynonymous SNV")){
      # if (variant_type == "LoF"){
      #   target.ix <- which(CNV_SNV_table$event_freq < freq &
      #                        CNV_SNV_table$effect_priority =="LoF" &
      #                        CNV_SNV_table$entrez_id %in% gsMain[[gsi]])
      # }
      # if (variant_type == "nonsynonymous"){
      #   target.ix <- which(CNV_SNV_table$event_freq < freq &
      #                        CNV_SNV_table$effect_priority =="nonsynonymous SNV" &
      #                        CNV_SNV_table$entrez_id %in% gsMain[[gsi]])
      # }
      target.ix <- which(CNV_SNV_table$event_freq < freq &
                           CNV_SNV_table$effect_priority == variant_type &
                           CNV_SNV_table$entrez_id %in% gsMain[[gsi]])
      
      target.child <- sum(CNV_SNV_table$Relation[target.ix] == "proband or affected sibling")
      target.parent <- sum(CNV_SNV_table$Relation[target.ix] == "parent")
      
      bg.child <- sum(CNV_SNV_table$Relation[-target.ix] == "proband or affected sibling") #inverted set of what is in target.ix
      bg.parent <- sum(CNV_SNV_table$Relation[-target.ix] == "parent")
      
      df.test <- data.frame("target" = c(target.child, target.parent-target.child),
                            "background" = c(bg.child, bg.parent-bg.child))
      
      test <- fisher.test(df.test, alternative = "greater")
      test.out <- rbind(test.out, 
                        data.frame("gene_set" = names(gsMain[gsi]),
                                   "event_freq" = format(freq, scientific = F),
                                   "variant_type" = variant_type,
                                   "OR" = signif(test$estimate, digits = 3),
                                   "P" = signif(test$p.value, digits = 3),
                                   target.child,
                                   target.parent,
                                   bg.child,
                                   bg.parent))
    }
  }
}

write.table(test.out, "MSSNG_gene_enrichment_analysis.tsv", sep="\t", row.names=F, quote=F, col.names=T)


