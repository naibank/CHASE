library(data.table)


## import pRecs
lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

## import metadata
metadata <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", data.table = F)
metadata$Relation <- gsub("other sibling", "unaffected sibling", metadata$Relation)
filtered_metadata <- read.table("/Users/shaniawu/SickKids CHASE/SSC/data/2021.11.02_SSC_filtered_metadata.txt")
metadata <- metadata[which(metadata$`Sample ID` %in% filtered_metadata$V1),] 

## find which IDs are parental, proband, or unaffected sibling
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation == "proband")]
meta_unaffectedsibID <- metadata$`Sample ID`[which(metadata$Relation == "unaffected sibling")]


####################################
########## Proband Table ###########
####################################

proband_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
CNV_SNV_table <- proband_table

## add relation column 
CNV_SNV_table$Relation <- NA
CNV_SNV_table$Relation[which(CNV_SNV_table$`#Sample` %in% meta_parentsID)] <- "parent"
CNV_SNV_table$Relation[which(CNV_SNV_table$`#Sample` %in% meta_probandID)] <- "proband"
CNV_SNV_table <- CNV_SNV_table[which(!CNV_SNV_table$Relation == "NA"),]

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
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"),])
        }
      }
      ## LOF SNV count
      if (count_table[row,"variant"] == "LoF"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband"&
                                                          table$effect_priority == "LoF"),])###Bank
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "LoF"),])###Bank
        }
      }
      ## Synonymous SNV count
      if (count_table[row,"variant"] == "Synonymous SNVs"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband"&
                                                          table$effect_priority == "synonymous SNV"),])
        }
        if (count_table[row,"relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "synonymous SNV"),])
        }
      }
      ## Nonsynonymous SNV count
      if (count_table[row,"variant"] == "Nonsynonymous SNVs"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband"&
                                                          table$effect_priority == "nonsynonymous SNV"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "nonsynonymous SNV"),])
        }
      }
      ## damaging missense SNV count
      if (count_table[row,"variant"] == "Damaging missenses"){
        if (count_table[row, "relation"] == "Probands"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband"&
                                                          table$damaging_missense_count > 3),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$damaging_missense_count > 3),])
        }
      }
    }
    
    plot_name =  paste("SSC_proband_eventfreq", freq, "_pRec", pRecx, sep="")###Bank
    
    ## plot count_table
    plot <- ggplot(count_table, aes(x = variant, y = count, fill = relation)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_label(aes(label = count), position = position_dodge(width = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(plot_name)  ###Bank
    
    ggsave(sprintf("%s.png", plot_name), width = 10)
  }
}




####################################
#### Unaffected Sibling Table ######
####################################

unaffectedsib_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
CNV_SNV_table <- unaffectedsib_table

## add relation column 
CNV_SNV_table$Relation <- NA
CNV_SNV_table$Relation[which(CNV_SNV_table$`#Sample` %in% meta_parentsID)] <- "parent"
CNV_SNV_table$Relation[which(CNV_SNV_table$`#Sample` %in% meta_unaffectedsibID)] <- "unaffected sibling"
CNV_SNV_table <- CNV_SNV_table[which(!CNV_SNV_table$Relation == "NA"),]

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
          count_table[row, "count"] <- nrow(table[which(table$Relation == "unaffected sibling"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"),])
        }
      }
      ## LOF SNV count
      if (count_table[row,"variant"] == "LoF"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "unaffected sibling"&
                                                          table$effect_priority == "LoF"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "LoF"),])
        }
      }
      ## Synonymous SNV count
      if (count_table[row,"variant"] == "Synonymous SNVs"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "unaffected sibling"&
                                                          table$effect_priority == "synonymous SNV"),])
        }
        if (count_table[row,"relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "synonymous SNV"),])
        }
      }
      ## Nonsynonymous SNV count
      if (count_table[row,"variant"] == "Nonsynonymous SNVs"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "unaffected sibling"&
                                                          table$effect_priority == "nonsynonymous SNV"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$effect_priority == "nonsynonymous SNV"),])
        }
      }
      ## damaging missense SNV count
      if (count_table[row,"variant"] == "Damaging missenses"){
        if (count_table[row, "relation"] == "Unaffected Sibling"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "unaffected sibling"&
                                                          table$damaging_missense_count > 3),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                          table$damaging_missense_count > 3),])
        }
      }
    }
    
    ## plot count_table
    plot_name =  paste("SSC_unaffectedsibling_eventfreq", freq, "_pRec", pRecx, sep="")###Bank
    
    plot <- ggplot(count_table, aes(x = variant, y = count, fill = relation)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_label(aes(label = count), position = position_dodge(width = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(plot_name)###Bank
    
    ggsave(sprintf("%s.png", plot_name), width = 10)
  }
}
