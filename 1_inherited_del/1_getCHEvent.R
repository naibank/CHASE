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
library(R.utils)
# library(writexl)
# library(openxlsx)
source("functions.R")

args = commandArgs(trailingOnly=TRUE)

print("read CNVs")  
datasets <- c("MSSNG_CG", "MSSNG_ILMN", paste("SPARK_WGS_", 1:4, sep=""), "SSC")
i <- which(datasets == args[1])

cnv_files <- list.files("../prerun_family_QCs/", full.names = T,  pattern = "mosaic.tagged.tsv")

snv_files <- c("MSSNG/CG/variants/SNVs+indels/exonic+splicing/",
               "MSSNG/ILMN/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_1/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_2/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_3/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_4/variants/SNVs+indels/exonic+splicing/",
               "SSC/variants/SNVs+indels/exonic+splicing/")


meta_files <- list.files("../prerun_family_QCs/", full.names = T, pattern = "metadata.tsv")
metatable <- read.delim(meta_files[i], stringsAsFactors = F)

# Getting CNVs (with different MAF limits) for samples sequenced by ILMN (only CNVs_10Percent is used further on)
CNVs_1Percent <- CNVfilter.procedure(CNVfile.path = cnv_files[i], metatable)


print(sprintf("getting CH events proband %s", datasets[i]))
Get_CH_Events(probands = metatable[which((!metatable$Relation %in% c("father", "mother")) & 
                                           metatable$Affection == 2 & metatable$family_member == "trios"), ],
             CNVdf = CNVs_1Percent,
             SNVfolder = snv_files[i]) %>%
  write_yaml(., sprintf("%s_Proband_CH_Data_CNV10P_SNV.yaml", datasets[i]))

print(sprintf("getting CH events unaff %s", datasets[i]))
Get_CH_Events(probands = metatable[which((!metatable$Relation %in% c("father", "mother")) & 
                                           metatable$Affection == 1 & metatable$family_member == "trios"), ],
              CNVdf = CNVs_1Percent,
              SNVfolder = snv_files[i]) %>%
  write_yaml(., sprintf("%s_Unaff_CH_Data_CNV10P_SNV.yaml", datasets[i]))
