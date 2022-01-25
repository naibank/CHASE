library(data.table)

cnv <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.CNVs.freq.10perc.tsv", data.table = F)

## Get list of CNVs that impact CDS:
cnv <- cnv[which(cnv$cds_symbol != ""),]


## Get list of all SNVs:
snv_1 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_1.tsv", data.table = F)
snv_2 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_2.tsv", data.table = F)
snv_3 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_3.tsv", data.table = F)
snv_4 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_4.tsv", data.table = F)
snv_5 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_5.tsv", data.table = F)
snv <- rbind(rbind(rbind(rbind(snv_1,snv_2),snv_3),snv_4),snv_5)
rm(list = paste("snv", 1:5, sep="_"))

## Get filtered metadata file:
metadata <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG_metadata.tsv")
filtered_metadata <- read.table("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/SickKids CHASE/2021.10.25_filtered_metadata_BE.txt")
metadata <- metadata[which(metadata$`Sample ID` %in% filtered_metadata$V1),] 
#list of samples belonging to families with a proband and both parents with SNV data

## get the list of genes and filter SNVs/Indels data based on the filtered list
cnv.genes <- unique(strsplit(paste(cnv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
snv <- snv[which(snv$gene_symbol %in% cnv.genes), ]
## Make tables:

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

write.table(table,"CNV_SNV_table.tsv", sep="\t", row.names=F, quote=F, col.names=T)

