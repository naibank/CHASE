#Required Packages
#BiocManager::install("GenomicRanges")
#BiocManager::install("Repitools")

library(GenomicRanges)
library(yaml)
library(dplyr)

options(stringsAsFactors = F)

source("functions.R")

datasets <- c("MSSNG_CG", "MSSNG_ILMN", paste("SPARK_WGS_", 1:4, sep=""), "SSC")

if(!dir.exists("SNV_in_CH_tables")){
  dir.create("SNV_in_CH_tables")
}

for(kid in c("Proband", "Unaff")){
  dt.out <- data.frame()
  for(dataset in datasets){
    message(dataset)
    dt <- Get_SNVs(sprintf("%s_%s_CH_Data_CNV10P_SNV.yaml", dataset, kid))
    if(nrow(dt.out) == 0)
      dt.out <- dt
    else
      dt.out <- rbind(dt.out[, intersect(names(dt.out), names(dt))], 
                      dt[, intersect(names(dt.out), names(dt))])
  }
  
  write.table(dt.out, sprintf("SNV_in_CH_tables/%s_parents.tsv", kid), sep="\t", row.names=F, quote=F, col.names=T)
}

  
  
