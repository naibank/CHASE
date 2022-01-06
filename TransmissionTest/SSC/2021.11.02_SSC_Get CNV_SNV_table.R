library(data.table)

cnv.main <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC.CNVs.freq.10perc.tsv", data.table = F)

## Get list of CNVs that impact CDS:
cnv.main <- cnv.main[which(cnv.main$cds_symbol != ""),]


## Get list of all SNVs:
snv_1 <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC.SNVs.freq.10perc_1.tsv", data.table = F)
snv_2 <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC.SNVs.freq.10perc_2.tsv", data.table = F)
snv_3 <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC.SNVs.freq.10perc_3.tsv", data.table = F)
snv_4 <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC.SNVs.freq.10perc_4.tsv", data.table = F)
snv_5 <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC.SNVs.freq.10perc_5.tsv", data.table = F)
snv <- rbind(rbind(rbind(rbind(snv_1,snv_2),snv_3),snv_4),snv_5)
rm(list = paste("snv", 1:5, sep="_")) 


## Get filtered metadata file:
metadata <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", data.table = F)
metadata$Relation <- gsub("other sibling", "unaffected sibling", metadata$Relation)
filtered_metadata <- read.table("/Users/shaniawu/SickKids CHASE/SSC/data/2021.11.02_SSC_filtered_metadata.txt")
metadata <- metadata[which(metadata$`Sample ID` %in% filtered_metadata$V1)] 

for(group in c("proband", "unaffected sibling")){
  cnv <- cnv.main[cnv.main$sample %in% metadata$`Sample ID`[metadata$Relation == group],]

  ## get the list of genes and filter SNVs/Indels data based on the filtered list
  cnv.genes <- unique(strsplit(paste(cnv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
  snv <- snv[which(snv$gene_symbol %in% cnv.genes), ]
  
  ## Make table:
  table <- data.frame()
  
  for (cnv_row in 1:nrow(cnv)){
    message(sprintf("%s out of %s", cnv_row, nrow(cnv)))
    cnv.tmp <- cnv[cnv_row, names(cnv)[!names(cnv) %in% names(snv)]]
    inh <- cnv.tmp$Inheritance
    Lof_cds_symbol <- strsplit(cnv.tmp$cds_symbol, "\\|")[[1]] #list of cds_symbols split at |
    
    probandID <- cnv.tmp$sample
    
    ## get non-transmitting and transmitting parent ID
    if (inh == "Paternal") {
      nonInhParentID <- metadata$`Mother ID`[which(metadata$`Sample ID`==probandID)]
      InhParentID <- metadata$`Father ID`[which(metadata$`Sample ID`==probandID)]
    } else {  #i.e. if inh == "Maternal"
      nonInhParentID <- metadata$`Father ID`[which(metadata$`Sample ID`==probandID)]
      InhParentID <- metadata$`Mother ID`[which(metadata$`Sample ID`==probandID)]
    }
    
    ## get SNVs of probands, transmitting and non-transmitting parent that match gene
    snv.inh.parent <- snv[which(snv$`#Sample` == InhParentID &
                                  snv$gene_symbol %in% Lof_cds_symbol), ]
    snv.noninh.parent <- snv[which(snv$`#Sample` == nonInhParentID &
                                  snv$gene_symbol %in% Lof_cds_symbol), ]
    snv.proband <- snv[which(snv$`#Sample` == probandID &
                           snv$gene_symbol %in% Lof_cds_symbol), ]
    
    ## get SNPs that found in both parents, as well as homozygous variants in non-inh parent. Then, remove them from the analysis
    snps.to.remove <- c(intersect(snv.inh.parent$"#id", snv.noninh.parent$"#id"),
                        snv.noninh.parent$"#id"[which(snv.noninh.parent$OZYG %in% c("alt-alt", "hom-alt"))])
    snv.noninh.parent <- snv.noninh.parent[which(!snv.noninh.parent$"#id" %in% snps.to.remove), ]
    
    ### retain only SNVs/Indels found in inh parent
    snv.proband <- snv.proband[which(snv.proband$"#id" %in% snv.noninh.parent$"#id"), ]
    snv.tmp <- rbind(snv.noninh.parent, snv.proband)
    
    ## check if snv.tmp is empty
    if (nrow(snv.tmp) == 0) {
      next #move to next cnv row
    }
    snv.tmp[, names(cnv.tmp)] <- cnv.tmp
    table <- rbind(table, snv.tmp)
  }
  
  write.table(table,sprintf("SSC_%s_CNV_SNV_table.tsv", gsub(" ", "_", group)), sep="\t", row.names=F, quote=F, col.names=T)
}


## remove CRVs from tables
proband_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table.tsv")
unaffectedsib_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table.tsv")

ls1 <- fread("/Users/shaniawu/SickKids CHASE/2021.11.01/data/CRVs/MSSNG+SSC.ASD135_LoF.tsv")
ls2 <- fread("/Users/shaniawu/SickKids CHASE/2021.11.01/data/CRVs/MSSNG+SSC.CNVs.tsv")
ls_exclude <- unique(rbind(ls1, ls2))
metadata <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", data.table = F)
familyID_exclude <- metadata$`Family ID`[which(metadata$`Sample ID` %in% ls_exclude$Sample)]
sampleID_exclude <- metadata$`Sample ID`[which(metadata$`Family ID` %in% familyID_exclude)]

proband_table_noCRVs <- proband_table[which(!proband_table$`#Sample` %in% sampleID_exclude),]
unaffectedsib_table_noCRVs <- unaffectedsib_table[which(!unaffectedsib_table$`#Sample` %in% sampleID_exclude),]

write.table(proband_table_noCRVs, "SSC_proband_CNV_SNV_table_noCRV.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
write.table(unaffectedsib_table_noCRVs, "SSC_unaffected_sibling_CNV_SNV_table_noCRV.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
