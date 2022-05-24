library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/Volumes/T5 ExFAT/Co-op/S22 Co-op SickKids (Scherer Lab)/CHASE_stor/SPARK/SPARK_TDT/")

## Get SPARK metadata
SPARK_meta <- data.frame()

for (i in 1:3) {
  SPARK_metadata_path <- paste("./data/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  meta <- read.delim(SPARK_metadata_path, stringsAsFactors = F)
  meta$WGS_ver <- i # add col indicating SPARK WGS version number

  SPARK_meta <- rbind(SPARK_meta, meta)
}

## TODO:
# - only keep parents from CNV_SNV_table to be imported
# - only keep nonsyn SNVs
# - separate tables into transmitted and non-transmitted SNVs
# - import meta to get child sample ID (neta) from parent sample ID (imported table)


## TRANSMITTED SNVS ##
IDs_and_pos <- fread('./data/SPARK.TDT.parent-US.SNV1.CNV0.01.tsv', data.table=F)
IDs_and_pos <- IDs_and_pos[which(IDs_and_pos$Relation == "parent" & # remove probands
                                   IDs_and_pos$effect_priority == "nonsynonymous SNV"),] # only keep nonsyn snvs
IDs_and_pos <- IDs_and_pos[which(IDs_and_pos$transmitted),]
    
write_script <- function(text, path='./data/SPARK_TDT_igv_false_snvinsv_parent-US_transmitted_script.sh' , append=T){
  write(text, path, append=append)
}

for (i in 1:nrow(IDs_and_pos)) {
  write_script('new',)
  write_script('genome hg38')
  write_script('snapshotDirectory ./')
  
  row <- IDs_and_pos[i,]
  ID <- row$`#Sample`
  child_ID <- row$sample #IDs_and_pos$sample[which(IDs_and_pos$`#Sample` == ID)]
      # df$#Sample = parentID
      # df$sample = child ID 
  chr <- row$CHROM
  pos <- row$POS
  
  # Get WGS version number for X.Sample
  ver <- SPARK_meta$WGS_ver[which(SPARK_meta$Sample.ID == ID)]
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SPARK_WGS_%s/CRAM/%s.cram', ver, ID)
  write_script(sprintf('load %s',load_path))
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SPARK_WGS_%s/CRAM/%s.cram', ver, child_ID)
  write_script(sprintf('load %s',load_path))
  
  location <- sprintf('%s:%s-%s',chr,pos-100,pos+100)
  write_script(sprintf('goto %s',location))
  
  write_script('colorBy READ_STRAND')
  
  filename <- paste(ID,chr,pos,sep='_')
  write_script(sprintf('snapshot %s',paste(filename,'.png',sep='')))
}




## NON-TRANSMITTED SNVS ##
IDs_and_pos <- fread('./data/SPARK.TDT.parent-US.SNV1.CNV0.01.tsv', data.table=F)
IDs_and_pos <- IDs_and_pos[which(IDs_and_pos$Relation == "parent" & # remove probands
                                   IDs_and_pos$effect_priority == "nonsynonymous SNV"),] # only keep nonsyn snvs
IDs_and_pos <- IDs_and_pos[which(!IDs_and_pos$transmitted),]

write_script <- function(text, path='./data/SPARK_TDT_igv_false_snvinsv_parent-US_non-transmitted_script.sh', append=T){
  write(text, path, append=append)
}

for (i in 1:nrow(IDs_and_pos)) {
  write_script('new',)
  write_script('genome hg38')
  write_script('snapshotDirectory ./')
  
  row <- IDs_and_pos[i,]
  ID <- row$`#Sample`
  child_ID <- row$sample #IDs_and_pos$sample[which(IDs_and_pos$`#Sample` == ID)]
  # df$#Sample = parentID
  # df$sample = child ID 
  chr <- row$CHROM
  pos <- row$POS
  
  # Get WGS version number for X.Sample
  ver <- SPARK_meta$WGS_ver[which(SPARK_meta$Sample.ID == ID)]
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SPARK_WGS_%s/CRAM/%s.cram', ver, ID)
  write_script(sprintf('load %s',load_path))
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SPARK_WGS_%s/CRAM/%s.cram', ver, child_ID)
  write_script(sprintf('load %s',load_path))
  
  location <- sprintf('%s:%s-%s',chr,pos-100,pos+100)
  write_script(sprintf('goto %s',location))
  
  write_script('colorBy READ_STRAND')
  
  filename <- paste(ID,chr,pos,sep='_')
  write_script(sprintf('snapshot %s',paste(filename,'.png',sep='')))
}

