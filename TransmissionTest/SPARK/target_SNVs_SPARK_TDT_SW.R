library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/Volumes/T5 ExFAT/Co-op/S22 Co-op SickKids (Scherer Lab)/CHASE_stor/SPARK/SPARK_TDT/")

meta <- rbind(read.delim("./data/SPARK_WGS_1_metadata_relfixed.tsv", stringsAsFactors = F),
              read.delim("./data/SPARK_WGS_2_metadata_relfixed.tsv", stringsAsFactors = F),
              read.delim("./data/SPARK_WGS_3_metadata_relfixed.tsv", stringsAsFactors = F))


Get_Transmitted_Indication <- function(table, comparison){
  ### Get_Transmitted_Indication(table, comparison): returns table with added "transmitted" column 
  ###     indicating if SNV is transmitted (T) or non-transmitted (F)
  ### df, str -> df
  
  ## Add family.UID col
  # table <- merge(table, meta[, c("Family.ID", "Sample.ID")], by.x = "#Sample", by.y = "Sample.ID", all.x = T)
  
  table$family.UID <- paste(table$sample,
                            ".", table$`#id`, sep="")
  
  ## Get transmitted SNVs (present in proband)
  if (comparison == "parent-proband"){
    transmitted.snvs.child <- table[which(table$Relation == "proband or affected sibling"),]
    transmitted.snvs.fam <- unique(transmitted.snvs.child$family.UID)
    
    transmitted.snvs.parent <- table[which(table$Relation == "parent" &
                                             table$family.UID %in% transmitted.snvs.fam),]
                                            
  } else if (comparison == "parent-US"){
    transmitted.snvs.child <- table[which(table$Relation == "unaffected sibling"),]
    transmitted.snvs.fam <- unique(transmitted.snvs.child$family.UID)
    
    transmitted.snvs.parent <- table[which(table$Relation == "parent" &
                                             table$family.UID %in% transmitted.snvs.fam),]
  }
  #transmitted.snvs.UID <- c(transmitted.snvs.child$UID, transmitted.snvs.parent$UID)

  
  ## Add transmitted indication col
  #table$transmitted <- ifelse((table$UID %in% transmitted.snvs.UID), T, F)
  table$transmitted <- ifelse((table$family.UID %in% transmitted.snvs.fam), T, F)
  
  return(table)
}


## Make target SNV tables
for (comparison in c("parent-proband", "parent-US")){
  if (comparison == "parent-proband"){
    CNV_SNV_table <- fread("./data/SPARK_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
  }
  if (comparison == "parent-US"){
    CNV_SNV_table <- fread("./data/SPARK_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
  }

  #CNV_SNV_table <- CNV_SNV_table_full[CNV_SNV_table_full$freq_max <= snv_freq,]
  
  ## group "LoF" SNVs
  CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                        CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank
  
  lof <- read.delim("./data/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)
  CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > 0.9]),]
  CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$event_freq < 0.01),]
  
  ## add damaging missense indication col
  CNV_SNV_table$damaging_missense <- CNV_SNV_table$damaging_missense_count > 3
  
  ## add transmitted indication col
  CNV_SNV_table <- Get_Transmitted_Indication(CNV_SNV_table, comparison)
  
  write.table(CNV_SNV_table, sprintf("./data/SPARK.TDT.%s.SNV1.CNV0.01.tsv", comparison), sep="\t", row.names=F)
}


### Check SNV counts with plot counts
SNV.count.table <- data.frame()

for (comparison in c("parent-proband", "parent-US")){
  if (comparison == "parent-proband"){
    target.snv.table <- fread("./data/SPARK.TDT.parent-proband.SNV1.CNV0.01.tsv", data.table = F)
  }
  if (comparison == "parent-US"){
    target.snv.table <- fread("./data/SPARK.TDT.parent-US.SNV1.CNV0.01.tsv", data.table = F)
  }

  ## Make count.table.row for each variant type
  for (variant_type in c("damaging_missense", "LoF", "nonsynonymous SNV","synonymous SNV")){
    if (variant_type == "damaging_missense"){
      transmitted.count <- nrow(target.snv.table[which(target.snv.table$transmitted &
                                                         target.snv.table$damaging_missense),])
      non.transmitted.count <- nrow(target.snv.table[which(!target.snv.table$transmitted &
                                                             target.snv.table$damaging_missense),])
    }
    else{
      transmitted.count <- nrow(target.snv.table[which(target.snv.table$transmitted &
                                                         target.snv.table$effect_priority == variant_type),])
      non.transmitted.count <- nrow(target.snv.table[which(!target.snv.table$transmitted &
                                                             target.snv.table$effect_priority == variant_type),])
    }
    SNV.count.table.row <- data.frame("comparison" = comparison,
                                      "variant_type" = variant_type, 
                                      "transmitted.count" = transmitted.count,
                                      "non.transmitted.count" = non.transmitted.count)
    SNV.count.table <- rbind(SNV.count.table, SNV.count.table.row)
  }
}

print(SNV.count.table)


