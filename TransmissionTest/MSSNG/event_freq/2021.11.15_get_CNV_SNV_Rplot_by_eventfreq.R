library(data.table)
library(ggplot2)
CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table/CNV_SNV_table_eventfreq.tsv", data.table = F)
# CNV_SNV_table <- fread("../data/CNV_SNV_table_eventfreq.tsv", data.table = F)

## group LOF mutations
# for (row in 1:nrow(CNV_SNV_table)){
#   if (!CNV_SNV_table[row, "effect_priority"] %in% c("nonsynonymous SNV", "synonymous SNV")){
#     CNV_SNV_table[row, "effect_priority"] <- "LoF SNV"
#   }
# }
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                                                           CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank

## import pRecs
lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

## plot count by event_freq and pRec cut-offs

for(freq in c(0.005, 0.001, 0.0005, 0.0001)){ # event freq cut-offs
  table <- CNV_SNV_table[CNV_SNV_table$event_freq < freq, ]
  
  for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
    if (!pRecx == 0){ # no pRec cut-off for pRecx == 0
      table <- table[table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)], ]
    }

    ## make count table
    variant <- c("All", "All", "LoF", "LoF", "Synonymous SNVs","Synonymous SNVs",
                 "Nonsynonymous SNVs", "Nonsynonymous SNVs", "Damaging missenses","Damaging missenses") ###Bank
    relation <- c("Probands or affected siblings", "Parents")
    
    count_table <- data.table(variant, relation, count = 0)
    
    for (row in 1:nrow(count_table)){
      ## All count
      if (count_table[row,"variant"] == "All"){
        if (count_table[row, "relation"] == "Probands or affected siblings"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband or affected sibling"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"),])
        }
      }
      ## LOF SNV count
      if (count_table[row,"variant"] == "LoF"){###Bank
        if (count_table[row, "relation"] == "Probands or affected siblings"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband or affected sibling"&
                                                                table$effect_priority == "LoF"),])###Bank
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                                table$effect_priority == "LoF"),])###Bank
        }
      }
      ## Synonymous SNV count
      if (count_table[row,"variant"] == "Synonymous SNVs"){
        if (count_table[row, "relation"] == "Probands or affected siblings"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband or affected sibling"&
                                                                table$effect_priority == "synonymous SNV"),])
        }
        if (count_table[row,"relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                                table$effect_priority == "synonymous SNV"),])
        }
      }
      ## Nonsynonymous SNV count
      if (count_table[row,"variant"] == "Nonsynonymous SNVs"){
        if (count_table[row, "relation"] == "Probands or affected siblings"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband or affected sibling"&
                                                                table$effect_priority == "nonsynonymous SNV"),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                                table$effect_priority == "nonsynonymous SNV"),])
        }
      }
      ## damaging missense SNV count
      if (count_table[row,"variant"] == "Damaging missenses"){
        if (count_table[row, "relation"] == "Probands or affected siblings"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "proband or affected sibling"&
                                                                table$damaging_missense_count > 3),])
        }
        if (count_table[row, "relation"] == "Parents"){
          count_table[row, "count"] <- nrow(table[which(table$Relation == "parent"&
                                                                table$damaging_missense_count > 3),])
        }
      }
    }
    plot_name =  paste("MSSNG_eventfreq", freq, "_pRec", pRecx, sep="")###Bank
    
    ## plot count_table
    plot <- ggplot(count_table, aes(x = variant, y = count, fill = relation)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_label(aes(label = count), position = position_dodge(width = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  + ggtitle(plot_name)###Bank
    
    ggsave(sprintf("%s.png", plot_name), width = 10)
  }
}
