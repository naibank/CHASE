## Required packages
library(broom)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(stringr)
library(shiny)
library(scales)
library(tidyr)
library(readr)
library(yaml)
library(data.table)
library(base)
library(GenomicRanges)
# library(writexl)
# library(openxlsx)


## STEP 1: Locate rare CNV deletions from ASD probands that have inherited it or P_denovo and envelop exons 
# Loading in CNV databases + changing first column name to Sample.ID (necessary to keep identifying column name the same among all databases)

CNVfilter.procedure <- function(CNVfile.path){
  CNVs <- as_tibble(read.delim(CNVfile.path, stringsAsFactors = FALSE))
  
  # To have uniform column name across all databases for Sample ID
  colnames(CNVs)[1] <- "Sample.ID"
  
  # Selecting only CNV deletions
  CNVs <- CNVs %>% dplyr::filter(variantTypeAnn == "deletion")
  
  # Filtering out CNVs with "NA", ambiguous and one parent sequenced inheritance
  # CNVs <- CNVs  %>% dplyr::filter(!is.na(Inheritance)) 
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance != "Ambiguous")
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance != "One_parent_sequenced") 
  
  # Selecting only Paternally or Maternally Inherited CNVs
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance == "Paternal" | Inheritance == "Maternal")
  
  # Selecting CNVs which envelop at least one CDS region
  CNVs <- CNVs %>% dplyr::filter(cds_symbol != ""& !is.na(cds_symbol))
  
  # Selecting CNVs where the QC is "ok"
  CNVs <- CNVs %>% dplyr::filter(Sample_QC == "ok")
  
  # Filtering out CNVs located on the sex chromosomes
  CNVs <- CNVs %>%  dplyr::filter(CHROM != "chrX" & CHROM != "chrY")
  
  return(CNVs)
}

print("TEST 1")

CNVs_10Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/1000G/CNVs/CNVs.1000G.freq_1percent.HQR.tsv")

print("TEST 2")

print("TEST 3")

#################################################################################################################################################

Get_CH_Events <- function(CNVdf, SNVfolder){
  CH_Data <- list()
  files <- list.files(SNVfolder, full.names = T)
  for(i in 1:length(files)){
    tryCatch(
      expr = {
        snvs <- fread(files[i])[high_quality == T]
        if (nrow(snvs) > 0) {
          colnames(snvs)[colnames(snvs)=='Comment'][1] = 'Comment1'
          proband <- snvs$'#Sample'[[1]]
          probandCNVs <- CNVdf %>% dplyr::filter(Sample.ID == proband)
          
          
          CH_Data[[proband]]$CNVs <- probandCNVs
          
          cnv <- "CNVs"
          snvout <- "SNVs"
          
          if(nrow(CH_Data[[proband]][[cnv]]) > 0){
            cnvs.g <- GRanges(CH_Data[[proband]][[cnv]]$CHROM, IRanges(CH_Data[[proband]][[cnv]]$START, CH_Data[[proband]][[cnv]]$END), "*")
            snvs.g <- GRanges(snvs$CHROM, IRanges(snvs$POS, snvs$POS), "*")
            olap <- findOverlaps(snvs.g, cnvs.g)
            snvs <- snvs[unique(olap@from), ]
            snvs <- as_tibble(snvs) %>% mutate(LoF = ifelse(str_detect(effect_impact, "Stop Gain-High") == T |
                                                              str_detect(effect_impact, "Frameshift-High") == T | 
                                                              str_detect(effect_impact, "Splice Site_High") == T |
                                                              str_detect(typeseq_priority, 'exonic;splicing') == T, T, F))
            CH_Data[[proband]][[snvout]] <- snvs
          }else{
            CH_Data[[proband]][[snvout]] <- as_tibble(data.frame())
          }
        }

      },
      error = function(e){ 
        message(paste(c("Error reading file:", files[i]), collapse=" "))
        message(e)
      },
      warning = function(w){
        message(paste(c("Warning reading file:", files[i]), collapse=" "))
        message(w)
      },
      finally = {
        next
      }
    )
    
    message(sprintf("File no.: %s of %s, Called: %s", i, length(files), files[i]))
  }
  
  return(CH_Data)
}   

print("getting CNV probands CH data")
Get_CH_Events(CNVdf = CNVs_10Percent,
              SNVfolder = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/1000G/exonic_splicing") %>%
  write_yaml(.,"1000G_CH_Data_CNV1P_SNV.yaml")

print("Done")