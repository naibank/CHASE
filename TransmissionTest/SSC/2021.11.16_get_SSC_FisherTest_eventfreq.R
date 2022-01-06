library(data.table)
meta <- read.delim("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", stringsAsFactors = F)
meta$Relation <- gsub("other sibling", "unaffected sibling", meta$Relation)
meta$Relation <- gsub("father|mother", "parent", meta$Relation)

for (relation in c("proband", "unaffected_sibling")){
  if (relation == "proband"){
    CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)

    child = "proband" 
  }
  if (relation == "unaffected_sibling"){
    CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)

    child = "unaffected sibling"
  }
  #CNV_SNV_table <- merge(CNV_SNV_table, meta[, c("Sample.ID", "Relation")], by.x = "#Sample", by.y = "Sample.ID", all.x = T)
  ## group LOF variants in effect_priority
  CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                        CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank
  
  lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

  test.out <- data.frame()
  
  for(freq in c(0.01, 0.001, 0.005, 0.0005, 0.0001)){ # event freq cut-offs
    for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
      for(variant_type in c("damaging_missense", "LoF", "nonsynonymous","synonymous")){
        
        if (!pRecx == 0){ 
          if (variant_type == "damaging_missense"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                 CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$damaging_missense_count > 3)
          }
          if (variant_type == "LoF"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                 CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$effect_priority =="LoF")
          }
          if (variant_type == "nonsynonymous"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                 CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$effect_priority =="nonsynonymous SNV")
          }
          if (variant_type == "synonymous"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                 CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$effect_priority =="synonymous SNV")
          }
          
        }else{ # no pRec cut-off for pRecx == 0
          if (variant_type == "damaging_missense"){
            target.ix <- which(CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$damaging_missense_count > 3)
          }
          if (variant_type == "LoF"){
            target.ix <- which(CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$effect_priority =="LoF")
          }
          if (variant_type == "nonsynonymous"){
            target.ix <- which(CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$effect_priority =="nonsynonymous SNV")
          }
          if (variant_type == "synonymous"){
            target.ix <- which(CNV_SNV_table$event_freq < freq &
                                 CNV_SNV_table$effect_priority =="synonymous SNV")
          }
        }
        
        target.child <- sum(CNV_SNV_table$Relation[target.ix] == child)
        target.parent <- sum(CNV_SNV_table$Relation[target.ix] == "parent")
        
        bg.child <- sum(CNV_SNV_table$Relation[-target.ix] == child)
        bg.parent <- sum(CNV_SNV_table$Relation[-target.ix] == "parent")
        
        df.test <- data.frame("target" = c(target.child, target.parent-target.child),
                              "background" = c(bg.child, bg.parent-bg.child))
        
        test <- fisher.test(df.test, alternative = "greater")
        test.out <- rbind(test.out, 
                          data.frame("event_freq" = format(freq, scientific = F),
                                     "pRec" = pRecx,
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
  table_name <- paste("SSC_", relation, "_FisherTest_eventfreq.tsv", sep = "")
  write.table(test.out, table_name, sep="\t", row.names=F, quote=F, col.names=T)
}
  