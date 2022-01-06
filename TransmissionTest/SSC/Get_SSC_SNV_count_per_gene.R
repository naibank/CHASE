library(data.table)

CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_noCRV.tsv", data.table = F)

## add Relation column
metadata <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", data.table = F)
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation == "proband")]

CNV_SNV_table$Relation <- NA
CNV_SNV_table$Relation[which(CNV_SNV_table$`#Sample` %in% meta_parentsID)] <- "parent"
CNV_SNV_table$Relation[which(CNV_SNV_table$`#Sample` %in% meta_probandID)] <- "proband"

## exclude chrX SNVs
CNV_SNV_table <- CNV_SNV_table[which(!CNV_SNV_table$CHROM == "chrX"),]


## make SNV count table
lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

CNV_SNV_table <- CNV_SNV_table[CNV_SNV_table$freq_max < 0.01, ]
CNV_SNV_table  <- CNV_SNV_table[CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > 0.9], ]

SNV_genes <- unique(CNV_SNV_table$gene_symbol)

gene_table <- data.frame(gene = SNV_genes)
gene_table$parent_nonsyn <- NA
gene_table$parent_syn <- NA
gene_table$proband_nonsyn <- NA
gene_table$proband_syn <- NA

for (x in 1:nrow(gene_table)){
  gene.tmp <- gene_table[x,]
  SNVs_for_gene <-  CNV_SNV_table[which(CNV_SNV_table$gene_symbol == gene.tmp$gene),]
  
  gene_table[x,"parent_nonsyn"] <- nrow(SNVs_for_gene[which(SNVs_for_gene$Relation == "parent" &
                                                              SNVs_for_gene$effect_priority == "nonsynonymous SNV"),])
  gene_table[x,"parent_syn"] <- nrow(SNVs_for_gene[which(SNVs_for_gene$Relation == "parent" &
                                                           SNVs_for_gene$effect_priority == "synonymous SNV"),])
  gene_table[x,"proband_nonsyn"] <- nrow(SNVs_for_gene[which(SNVs_for_gene$Relation == "proband" &
                                                               SNVs_for_gene$effect_priority == "nonsynonymous SNV"),])
  gene_table[x,"proband_syn"] <- nrow(SNVs_for_gene[which(SNVs_for_gene$Relation == "proband" &
                                                            SNVs_for_gene$effect_priority == "synonymous SNV"),])
}

gene_table <- setorder(gene_table, -"proband_nonsyn")

write.table(gene_table, "SSC_SNV_count_per_gene.tsv", sep="\t", row.names=F, quote=F, col.names=T)

