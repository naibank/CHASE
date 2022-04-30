library(data.table)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

## import pRecs
lof <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

####################################
########## Proband Table ###########
####################################

proband_table <- fread("SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
CNV_SNV_table <- proband_table
#Adjust SNV freq
# CNV_SNV_table <- CNV_SNV_table[CNV_SNV_table$freq_max <= 0.01,]
## group LOF mutations
# for (row in 1:nrow(CNV_SNV_table)){
#   if (!CNV_SNV_table[row, "effect_priority"] %in% c("nonsynonymous SNV", "synonymous SNV")){
#     CNV_SNV_table[row, "effect_priority"] <- "LoF SNV"
#   }
# }
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                      CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF"###Bank

## plot count by event_freq and pRec cut-offs

for(freq in c(0.01, 0.005, 0.001, 0.0005, 0.0001)){ # event freq cut-offs
  table <- CNV_SNV_table[CNV_SNV_table$event_freq < freq, ]
  
  for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
    if (!pRecx == 0){ # no pRec cut-off for pRecx == 0
      table <- table[table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)], ]###Bank
    }
    
    ## make count table
    variant <- c("All", "All", "LoF", "LoF", "Synonymous SNVs","Synonymous SNVs",
                 "Nonsynonymous SNVs", "Nonsynonymous SNVs", "Damaging missenses","Damaging missenses")###Bank
    relation <- c("Probands", "Parents")
    
    count_table <- data.table(variant, relation, count = 0)
    
    for (row in 1:nrow(count_table)){
      ## All count
      if (count_table[row,"variant"] == "All"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "proband or affected sibling"), c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"), c("#Sample", "#id")]))
        }
      }
      ## LOF SNV count
      if (count_table[row,"variant"] == "LoF"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "proband or affected sibling"&
                                                          table$effect_priority == "LoF"),c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                                 table$effect_priority == "LoF"),c("#Sample", "#id")]))
        }
      }
      ## Synonymous SNV count
      if (count_table[row,"variant"] == "Synonymous SNVs"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "proband or affected sibling"&
                                                                 table$effect_priority == "synonymous SNV"),c("#Sample", "#id")]))
        }
        if (count_table[row,"relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent" &
                                                                 table$effect_priority == "synonymous SNV"),c("#Sample", "#id")]))
        }
      }
      ## Nonsynonymous SNV count
      if (count_table[row,"variant"] == "Nonsynonymous SNVs"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "proband or affected sibling"&
                                                                 table$effect_priority == "nonsynonymous SNV"),c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                                 table$effect_priority == "nonsynonymous SNV"),c("#Sample", "#id")]))
        }
      }
      ## damaging missense SNV count
      if (count_table[row,"variant"] == "Damaging missenses"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "proband or affected sibling"&
                                                          table$damaging_missense_count > 3),c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                          table$damaging_missense_count > 3),c("#Sample", "#id")]))
        }
      }
    }
    
    plot_name =  paste("../TT/Combined Data/SSC/SSC_proband_eventfreq", freq, "_pRec", pRecx, "_trate", sep="")###Bank
    
    trate_table = data.table(variant_type = unique(count_table$variant), trate = NA, 
                              n = count_table$count[which(count_table$relation == "Probands")])
    trate_table$trate <- count_table$count[which(count_table$relation == "Probands")] /
      count_table$count[which(count_table$relation == "Parents")]
    
    ## plot trate_table
    plot <- ggplot(trate_table, aes(x = variant_type, y = trate)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_label(aes(label = round(trate, digits = 2))) +
      geom_text(aes(label = paste("n = ", n), vjust = 3, color = "white"), show.legend = F) +
      labs(x = "Variant Type", y = "Transmisstion Rate") +
      geom_hline(yintercept = 0.5, lty = 2, color = "red") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.title = element_text(size = 12, face = "bold"))
    
    ggsave(sprintf("%s.png", plot_name), width = 5)
  }
}




####################################
#### Unaffected Sibling Table ######
####################################

unaffectedsib_table <- fread("SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
CNV_SNV_table <- unaffectedsib_table
# CNV_SNV_table <- CNV_SNV_table[CNV_SNV_table$freq_max <= 0.01,]

## group LOF mutations
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                      CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF"###Bank
# for (row in 1:nrow(CNV_SNV_table)){
#   if (!CNV_SNV_table[row, "effect_priority"] %in% c("nonsynonymous SNV", "synonymous SNV")){
#     CNV_SNV_table[row, "effect_priority"] <- "LoF SNV"
#   }
# }

## plot count by event_freq and pRec cut-offs

for(freq in c(0.01, 0.005, 0.001, 0.0005, 0.0001)){ # event freq cut-offs
  table <- CNV_SNV_table[CNV_SNV_table$event_freq < freq, ]
  
  for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
    if (!pRecx == 0){ # no pRec cut-off for pRecx == 0
      table <- table[table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)], ]
    }
    
    ## make count table
    variant <- c("All", "All", "LoF", "LoF", "Synonymous SNVs","Synonymous SNVs",
                 "Nonsynonymous SNVs", "Nonsynonymous SNVs", "Damaging missenses","Damaging missenses")###Bank
    relation <- c("Unaffected Sibling", "Parents")
    
    count_table <- data.table(variant, relation, count = 0)
    
    for (row in 1:nrow(count_table)){
      ## All count
      if (count_table[row,"variant"] == "All"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "unaffected sibling"), c("#Sample", "#id")])) 
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"), c("#Sample", "#id")]))
        }
      }
      ## LOF SNV count
      if (count_table[row,"variant"] == "LoF"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "unaffected sibling"&
                                                          table$effect_priority == "LoF"), c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "LoF"), c("#Sample", "#id")]))
        }
      }
      ## Synonymous SNV count
      if (count_table[row,"variant"] == "Synonymous SNVs"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "unaffected sibling"&
                                                          table$effect_priority == "synonymous SNV"),c("#Sample", "#id")]))
        }
        if (count_table[row,"relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "synonymous SNV"),c("#Sample", "#id")]))
        }
      }
      ## Nonsynonymous SNV count
      if (count_table[row,"variant"] == "Nonsynonymous SNVs"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "unaffected sibling"&
                                                          table$effect_priority == "nonsynonymous SNV"),c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "nonsynonymous SNV"),c("#Sample", "#id")]))
        }
      }
      ## damaging missense SNV count
      if (count_table[row,"variant"] == "Damaging missenses"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "unaffected sibling"&
                                                          table$damaging_missense_count > 3),c("#Sample", "#id")]))
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(unique(table[which(table$Relation == "parent"&
                                                          table$damaging_missense_count > 3),c("#Sample", "#id")]))
        }
      }
    }
    
    plot_name =  paste("../TT/Combined Data/SSC/SSC_unaffectedsibling_eventfreq", freq, "_pRec", pRecx, "_trate", sep="")
    
    trate_table <- data.table(variant_type = unique(count_table$variant), trate = NA, 
                              n = count_table$count[which(count_table$relation == "Unaffected Sibling")])
    trate_table$trate = count_table$count[which(count_table$relation == "Unaffected Sibling")] /
      count_table$count[which(count_table$relation == "Parents")]
    
    ## plot trate_table
    plot <- ggplot(trate_table, aes(x = variant_type, y = trate)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_label(aes(label = round(trate, digits = 2))) +
      geom_text(aes(label = paste("n = ", n), vjust = 3, color = "white"), show.legend = F) +
      labs(x = "Variant Type", y = "Transmisstion Rate") +
      geom_hline(yintercept = 0.5, lty = 2, color = "red") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.title = element_text(size = 12, face = "bold"))
    
    ggsave(sprintf("%s.png", plot_name), width = 5)
  }
}
