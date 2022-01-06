library(data.table)
library(ggplot2)
metadata <- fread("../data/MSSNG_metadata.tsv")

## find which IDs are parental and which are proband
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation == "proband")]
lof <- read.delim("data/gnomAD_pLI_v2.1_Feb8.txt", stringsAsFactors = F)
## separate table of probands and parents
for(freq in c(0.1, 0.05, 0.01, 0.001)){
  table <- fread("CNV_SNV_table.tsv")
  table <- table[table$freq_max < freq, ]
  table <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.9], ]
  ## if you like to add a column to indicate the relation, another way to do it is as below
  ## table$Relation <- NA
  ## table$Relation[which(table$`#Sample` %in% meta_parentsID) <- "parent"
  ## table$Relation[which(table$`#Sample` %in% meta_probandID) <- "proband"
  ## then, you can remove table with Relation is NA. Doing it this way will save a bit of memory, as you don't need parents_table and probands_table variables
  
  parents_table <- table[which(table$`#Sample` %in% meta_parentsID)]
  probands_table <- table[which(table$`#Sample` %in% meta_probandID)]
  
  ## add column indicating relation
  parents_table$Relation <- c("parent")
  probands_table$Relation <- c("proband")
  
  ## parents & proband table merged!
  table <- rbind(parents_table, probands_table)
  
  ## plot effect_priority vs. count
  plot <- ggplot(data = table, aes(x = effect_priority, fill = Relation))
  plot + geom_bar(stat = "count", position = position_dodge()) + geom_label(stat = "count", aes(label = ..count.., fill = Relation), 
                                                                           position = position_dodge(width = 1))
  
  ggsave(sprintf("%s.png", freq), width = 10)
}
