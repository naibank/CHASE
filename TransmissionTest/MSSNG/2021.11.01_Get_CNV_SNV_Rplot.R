library(data.table)
library(ggplot2)

metadata <- fread("MSSNG_metadata.tsv")

## find which IDs are parental and which are proband
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation %in% c("affected sibling", "proband"))]

meta_probandID1 <- metadata$`Sample ID`[which(metadata$Relation == "proband")]
meta_affectedsiblingID <- metadata$`Sample ID`[which(metadata$Relation == "affected sibling")]


lof <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

## import CNV_SNV_table with CRVs excluded
table <- fread("CNV_SNV_table.tsv")
table$V1 <- NULL
#table0.1 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.1], ]
#table0.01 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.01], ]
#table0.05 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.05], ]
#table0.001 <- table[table$gene_symbol %in% lof$gene[lof$pRec > 0.001], ]
#table0 <- table #[table$gene_symbol %in% lof$gene[lof$pRec > 0], ]

affectedsib_table <- table[which(table$`#Sample` %in% meta_affectedsiblingID),]
pro_table <- table[which(table$`#Sample` %in% meta_probandID1),]


for(freq in c(0.1, 0.05, 0.01, 0.001)){
  table <- fread("CNV_SNV_table.tsv")
  ## add a column to indicate the relation
  table$Relation <- NA
  table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
  table$Relation[which(table$`#Sample` %in% meta_probandID)] <- "proband or affected sibling"
  table <- table[which(!table$Relation == "NA"),]
  
  table <- table[which(!table$CHROM == "chrX")] #exclude chrX variants
  
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
    
    plot_name =  paste(freq, pRecx, sep="_sRec")
    
    ggsave(sprintf("../TT/MSSNG_%s.png", plot_name), width = 10)}
}


