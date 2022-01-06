library(data.table)
library(ggplot2)

metadata <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG_metadata.tsv")

## find which IDs are parental and which are proband
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation == "proband")]
lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

## separate table of probands and parents
for(freq in c(0.1, 0.05, 0.01, 0.001)){
  table <- fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table.tsv")
  table <- table[which(!table$CHROM == "chrX")]
  table <- table[table$freq_max < freq, ]
  
  for(pRecx in c(0, 0.7, 0.8, 0.9)){
    table <- table[table$gene_symbol %in% lof$gene[lof$pRec > pRecx], ]

    ## if you like to add a column to indicate the relation, another way to do it is as below
    #table$Relation <- NA
    #table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
    #table$Relation[which(table$`#Sample` %in% meta_probandID)] <- "proband"
    ## then, you can remove table with Relation is NA. Doing it this way will save a bit of memory, as you don't need parents_table and probands_table variables
    #table <- table[which(table$Relation == "NA"),]
  
    parents_table <- table[which(table$`#Sample` %in% meta_parentsID)]
    probands_table <- table[which(table$`#Sample` %in% meta_probandID)]
    ### add column indicating relation
    parents_table$Relation <- c("parent")
    probands_table$Relation <- c("proband")
    ### parents & proband table merged!
    table <- rbind(parents_table, probands_table)
  
    ## plot effect_priority vs. count
    plot <- ggplot(data = table, aes(x = effect_priority, fill = Relation))
    plot + geom_bar(stat = "count", position = position_dodge()) + 
      geom_label(stat = "count", aes(label = ..count.., fill = Relation), position = position_dodge(width = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_name =  paste(freq, pRecx, sep="_sRec")
    
    ggsave(sprintf("%s.png", plot_name), width = 10)}
}

