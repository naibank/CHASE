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

setwd("/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MGRB/data_processing/data/")

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


CNVs_10Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MGRB/CNVs/CNVs.MGRB.freq_1percent.HQR.tsv")


#################################################################################################################################################

Get_CH_Events <- function(CNVdf, SNVfolder){
  CH_Data <- list()
  files <- list.files(SNVfolder, full.names = T)
  for(i in 1:length(files)){
    snvs <- fread(files[i])
    message(i)
    tryCatch(
      expr = {
        snvs <- snvs[high_quality == T]
        if (nrow(snvs) > 0) {
          colnames(snvs)[colnames(snvs)=='Comment'][1] = 'Comment1'
          
          ## add "#Sample" column
          file.name <- basename(files[i])
          sample.ID <- substring(file.name, 1, 5)
          snvs$"#Sample" <- sample.ID
          snvs <- snvs[,c(216, 1:215)] # move "#Sample" col to front
          
          ## change "CHROM" col format (e.g. "1" -> "chr1")
          snvs$`#CHROM` <- paste("chr", snvs$`#CHROM`, sep="")
          
          ## remove sample ID prefixes from columns 20-30
          names(snvs)[20:30] <- substring(names(snvs)[20:30], 7)
          
          ## remove extra "DP" column
          snvs <- snvs[,-11]
          
          ## add "alt_fraction" column to snvs
          snvs$"alt_fraction" <- as.numeric(snvs$AD_ALT) / as.numeric(snvs$DP) # alt_fraction = AD_ALT / DP
          
          proband <- snvs$'#Sample'[[1]]
          probandCNVs <- CNVdf %>% dplyr::filter(Sample.ID == proband)
          
          CH_Data[[proband]]$CNVs <- probandCNVs
          
          cnv <- "CNVs"
          snvout <- "SNVs"
          
          if(nrow(CH_Data[[proband]][[cnv]]) > 0){
            cnvs.g <- GRanges(CH_Data[[proband]][[cnv]]$CHROM, IRanges(CH_Data[[proband]][[cnv]]$START, CH_Data[[proband]][[cnv]]$END), "*")
            snvs.g <- GRanges(snvs$`#CHROM`, IRanges(snvs$start, snvs$end), "*")
            olap <- findOverlaps(snvs.g, cnvs.g)
            snvs <- snvs[unique(olap@from), ]
            snvs <- as_tibble(snvs) %>% mutate(LoF = ifelse(str_detect(effect_priority, "frameshift deletion") == T |
                                                              str_detect(effect_priority, "frameshift intersetion") == T |
                                                              str_detect(effect_priority, "stopgain") == T |
                                                              str_detect(typeseq_priority, "splicing") == T |
                                                              str_detect(typeseq_priority, "exonic;splicing") == T , T, F))
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


Get_CH_Events(CNVdf = CNVs_10Percent,
              SNVfolder = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MGRB/exonic+splicing") %>%
  write_yaml(.,"MGRB_CH_Data_CNV1P_SNV.yaml")

