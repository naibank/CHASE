library(data.table)
library(ggplot2)

metadata <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", data.table = F)
metadata$Relation <- gsub("other sibling", "unaffected sibling", metadata$Relation)
filtered_metadata <- read.table("/Users/shaniawu/SickKids CHASE/SSC/data/2021.11.02_SSC_filtered_metadata.txt")
metadata <- metadata[which(metadata$`Sample ID` %in% filtered_metadata$V1),] 

## find which IDs are parental and which are proband
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation == "proband")]
meta_unaffectedsib <- metadata$`Sample ID`[which(metadata$Relation == "unaffected sibling")]

table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_CNV_SNV_table.tsv", data.table = F)
#table0.1 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.1], ]
#table0.01 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.01], ]
#table0.05 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.05], ]
#table0.001 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.001], ]
#table0 <- table #[table$gene_symbol %in% lof$gene[lof$pRec > 0], ]

lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

proband_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table.tsv", data.table = F)
unaffectedsib_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table.tsv", data.table = F)

## plot proband table
for(freq in c(0.1, 0.05, 0.01, 0.001)){
  table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_noCRV.tsv", data.table = F)
  ## add a column to indicate the relation
  table$Relation <- NA
  table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
  table$Relation[which(table$`#Sample` %in% meta_probandID)] <- "proband"
  table <- table[which(!table$Relation == "NA"),]
  
  table <- table[which(!table$CHROM == "chrX"),] #exclude chrX variants
  table <- table[table$freq_max < freq, ]
  
  for(pRecx in c(0, 0.7, 0.9)){
    if (!pRecx == 0){
      table <- table[table$gene_symbol %in% lof$gene[lof$pRec > pRecx], ]
      }
    message(sprintf("freq%fpRec%f", freq, pRecx))
    ## plot effect_priority vs. count
    plot <- ggplot(data = table, aes(x = effect_priority, fill = Relation))
    plot + geom_bar(stat = "count", position = position_dodge()) + 
      geom_label(stat = "count", aes(label = ..count.., fill = Relation), position = position_dodge(width = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_name =  paste("SSC_proband_", freq, paste("_sRec",pRecx, sep = ""), "_noCRV", sep = "")
    
    ggsave(sprintf("%s.png", plot_name), width = 10)}}


## plot unaffected sibling table
for(freq in c(0.1, 0.05, 0.01, 0.001)){
  table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table_noCRV.tsv", data.table = F)
  ## add a column to indicate the relation
  table$Relation <- NA
  table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
  table$Relation[which(table$`#Sample` %in% meta_unaffectedsib)] <- "unaffected sibling"
  table <- table[which(!table$Relation == "NA"),]
  
  table <- table[which(!table$CHROM == "chrX"),] #exclude chrX variants
  table <- table[table$freq_max < freq, ]
  
  for(pRecx in c(0, 0.7, 0.9)){
    if (!pRecx == 0){
      table <- table[table$gene_symbol %in% lof$gene[lof$pRec > pRecx], ]
    }
    message(sprintf("freq%fpRec%f", freq, pRecx))
    ## plot effect_priority vs. count
    plot <- ggplot(data = table, aes(x = effect_priority, fill = Relation))
    plot + geom_bar(stat = "count", position = position_dodge()) + 
      geom_label(stat = "count", aes(label = ..count.., fill = Relation), position = position_dodge(width = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_name =  paste("SSC_unaffectedsibling_", freq, paste("_sRec",pRecx, sep = ""), "_noCRV", sep = "")
    
    ggsave(sprintf("%s.png", plot_name), width = 10)}}


