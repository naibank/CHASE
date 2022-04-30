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
library(GenomicRanges)

#setwd("/hpf/largeprojects/tcagstor/tcagstor_tmp/dmager/SSC")
# setwd("/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/SSC")

## STEP 1: Locate rare CNV deletions of ASD probands, who have inherited it or P_denovo and envelop exons 
CNVfilter.procedure <- function(CNVfile.path){
  CNVs <- as_tibble(read.delim(CNVfile.path,stringsAsFactors = F))
  # To have uniform column name across all databases for Sample ID
  colnames(CNVs)[1] <- "Sample.ID"
  
  # Selecting only CNV deletions
  CNVs <- CNVs %>% dplyr::filter(variantTypeAnn == "deletion")
  
  # Filtering out CNVs with "NA", Ambiguous and one parent sequenced inheritance
  # CNVs <-CNVs  %>% dplyr::filter(!is.na(Inheritance)) 
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance != "Ambiguous")
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance != "One_parent_sequenced") 
  # 
  # Selecting only Paternally or Maternally Inherited CNVs
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance == "Paternal"|Inheritance == "Maternal")
  
  # Selecting CNVs which envelop at least one CDS region
  CNVs <- CNVs %>% dplyr::filter(cds_symbol != ""& !is.na(cds_symbol))
  
  # Selecting CNVs where the QC is "ok"
  CNVs <- CNVs %>% dplyr::filter(Sample_QC == "ok")
  
  # Filtering out CNVs located on the sex chromosomes
  CNVs <- CNVs %>%  dplyr::filter(CHROM != "chrX" & CHROM != "chrY")
  
  return(CNVs)
}

print("TEST 1")

# Loading in CNV databases + changing first column name to Sample.ID
# Applying CNV filter procedure to SSC CNV 10P database
CNVs_10Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/CNVs/CNVs.SSC.freq_10percent.HQR.IntFreq.tsv")

## STEP 2
# Family to exclude because proband also has an affected mother
Family_to_exclude <- as_tibble(read.delim("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/SSC_metadata.tsv", stringsAsFactors = F)) %>%
  filter(Relation == "proband") %>%
  filter(Sample.ID %in% Mother.ID | Sample.ID %in% Father.ID) %>% pull(Family.ID)

# Subjects to exclude that failed CNV QC (Subjects can be added if needed)
Subjects_to_exclude <- c("SSC09651","SSC05770","SSC09923","SSC04868","SSC10162","SSC07001","SSC07819","SSC09280","SSC08169","SSC10798","SSC11555","SSC02241",
                         "SSC05551","SSC08065","SSC11881","SSC10250","SSC05055","SSC11874","SSC05183","SSC05809","SSC07802","SSC01113","SSC09431","SSC07287",
                         "SSC10874","SSC03297","SSC06244","SSC12781","SSC07579","SSC09671","SSC00091","SSC02016","SSC09360","SSC06759","SSC00089","SSC12226",
                         "SSC03785","SSC11924","SSC10173","SSC00315","SSC07530","SSC05498","SSC01947","SSC11590","SSC07000","SSC02716","SSC11250","SSC07806",
                         "SSC12399","SSC11356","SSC00671","SSC03166","SSC02124","SSC07797","SSC12407","SSC02666","SSC02033","SSC10881","SSC12196","SSC06563",
                         "SSC09092","SSC05299","SSC01117","SSC06572","SSC11606","SSC09929","SSC04958","SSC02933","SSC02929","SSC07759","SSC01121","SSC04475",
                         "SSC12763","SSC10621","SSC09472","SSC10759","SSC06390","SSC06533","SSC12459","SSC02132","SSC01128","SSC03304","SSC08188","SSC05990",
                         "SSC08521","SSC11530","SSC11712","SSC01120","SSC07675","SSC10098","SSC06580","SSC07761","SSC09098","SSC07131","SSC02188","SSC06103",
                         "SSC05559","SSC08689","SSC04278")


# Locate probands ID's
LocateProbands <- function(metadata.filepath,CNV.database,Family_to_exclude,Subjects_to_exclude){
  probands <- as_tibble(read.delim(metadata.filepath, stringsAsFactors = FALSE)) %>%
    dplyr::filter(Relation == "proband") %>%
    filter(Exclude.because.re.sequenced. != "yes") %>%
    filter(!Family.ID %in% Family_to_exclude &
             !Sample.ID %in% Subjects_to_exclude) #Exclude probands if there sample.ID is located in the vector of individuals that need to be excluded
  
  probands <- probands %>%
    dplyr::filter(Sample.ID %in% CNV.database$Sample.ID)
  return(probands)
}

Locate_Usiblings <- function(metadata.filepath, CNV.database, Family_to_exclude, Subjects_to_exclude){
  probands <- as_tibble(read.delim(metadata.filepath, stringsAsFactors = FALSE)) %>%
    dplyr::filter(Relation %in% c("unaffected sibling", "other sibling")) %>%
    filter(Exclude.because.re.sequenced. != "yes") %>%
    filter(!Family.ID %in% Family_to_exclude &
             !Sample.ID %in% Subjects_to_exclude) #Exclude probands if there sample.ID is located in the vector of individuals that need to be excluded
  
  probands <- probands %>%
    dplyr::filter(Sample.ID %in% CNV.database$Sample.ID)
  return(probands)
}

# To get probands who inherited at least one CNV deletion from a parent. 
probands_10P <- LocateProbands(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/SSC_metadata.tsv",
                               CNV.database = CNVs_10Percent,
                               Family_to_exclude = Family_to_exclude,
                               Subjects_to_exclude = Subjects_to_exclude)

# To get the unaffected siblings who Inherited at least one CNV deletion from a parent. 
Unaffected_siblings10P <- Locate_Usiblings(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/SSC_metadata.tsv",
                                         CNV.database = CNVs_10Percent,
                                         Family_to_exclude = Family_to_exclude,
                                         Subjects_to_exclude = Subjects_to_exclude)

# Locate parents of those probands
LocateParents <- function(metadata.filepath, probands.data, Family_to_exclude, Subjects_to_exclude){
  # Filters out parents/samples that have been resequenced, are not in the vectors of samples that need to be excluded and keeps only the parents of the probands with at least one inherited CNV deletion.
  parents <- as_tibble(read.delim(metadata.filepath, stringsAsFactors = FALSE)) %>%
    filter(Exclude.because.re.sequenced. != "yes") %>%
    filter(!Family.ID %in% Family_to_exclude &
             !Sample.ID %in% Subjects_to_exclude) %>%
    dplyr::filter(Sample.ID %in% probands.data$Father.ID | Sample.ID %in% probands.data$Mother.ID)
  return(parents)
}

parents_10P <- LocateParents(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/SSC_metadata.tsv",
                             probands.data = probands_10P,
                             Family_to_exclude = Family_to_exclude,
                             Subjects_to_exclude = Subjects_to_exclude)

print("TEST 2")

## STEP 3: A loop that looks per proband if there are any CH events, gets relevant CNV + SNV data (from proband, mother and father) if so and figures out the inheritance of SNV
# To store data

# Function and loop to find possible CH events in the CNV databases

Get_CH_Events <- function(probands, CNVdf, SNVfolder){
  CH_Data <- list()
  
  for(i in 1:nrow(probands)){
    
    #Get the filepath to the proband's relevant SNV file
    probandFilepath <- str_c(c("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/", SNVfolder, "/", 
                               probands$Sample.ID[i], ".tsv.gz"), collapse = "") # Search proband SNV data
    
    proband <- probands$Sample.ID[i]
    
    # Empty vector for probands' mother
    mother <- vector("character", length = 1)
    
    # Get the filepath to mother SNV file
    mother <- probands %>% dplyr::filter(probands$Sample.ID == probands$Sample.ID[i]) %>% 
      pull(Mother.ID) #search probands' mother SNV data
    motherFilepath <- str_c(c("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/", SNVfolder, "/",
                              mother[1], ".tsv.gz"), collapse = "")
    
    # Empty vector for probands' father
    father <- vector("character", length = 1)
    
    # Get the filepath to father SNV file
    father <- probands %>% dplyr::filter(probands$Sample.ID == probands$Sample.ID[i]) %>%
      pull(Father.ID) #search probands' father SNV data
    fatherFilepath <- str_c(c("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SSC/", SNVfolder, "/",
                              father[1], ".tsv.gz"), collapse = "")
    
    # Skip probands with no parent ID information
    if (probands$Sample.ID == '-' | mother[1] == '-' | father[1] == '-') {
      next
    }
    
    # Identify if proband has a CH event
    # METHOD 1: Identifies all SNVs that are within a CNV region of a proband (1) --> therefore, SNVs are on the allele without the deletion (CNV)
    # METHOD 2: Identifies all SNVs that are somewhat outside of CNV boundaries (2) --> however, located on the other allele + in a gene affected by the CNV
    
    # Load CNV data from proband, mother and father
    probandCNVs <- CNVdf %>% dplyr::filter(Sample.ID == proband)
    motherCNVs <- CNVdf %>% dplyr::filter(Sample.ID == mother)
    fatherCNVs <- CNVdf %>% dplyr::filter(Sample.ID == father)
    
    # get CNVs data per proband (mother and father CNVs as well)
    # if inheritance CNVs is not found in parent, we removed
    # filter out non transmitted CNVs
    # probandFilteredCNVs <- data.frame()
    # for(inheritance in c("Maternal", "Paternal")){
    #   parentCNVs <- motherCNVs
    #   parent <- "Mother"
    #   if(inheritance == "Paternal"){
    #     parentCNVs <- fatherCNVs
    #     parent <- "Father"
    #   }
    #   if(sum(probandCNVs$Inheritance == inheritance) > 0){
    #     inheritanceCNVs <- probandCNVs %>% dplyr::filter(Inheritance == inheritance)
    #     
    #     
    #     inheritanceCNVs.g <- GRanges(inheritanceCNVs$CHROM, IRanges(inheritanceCNVs$START, inheritanceCNVs$END), "*")
    #     if(nrow(parentCNVs) > 0){
    #       parentCNVs.g <- GRanges(parentCNVs$CHROM, IRanges(parentCNVs$START, parentCNVs$END), "*")
    #       olap <- findOverlaps(inheritanceCNVs.g, parentCNVs.g)
    #       parentCNVs <- parentCNVs[unique(olap@to), ]
    #       probandFilteredCNVs <- rbind(probandFilteredCNVs, inheritanceCNVs[unique(olap@from), ])
    #     }
    #     
    #     CH_Data[[proband]][[sprintf("%s.CNV", parent)]] <- parentCNVs
    #   }else{
    #     CH_Data[[proband]][[sprintf("%s.CNV", parent)]] <- as_tibble(data.frame())
    #   }
    # }
    # CH_Data[[proband]]$CNVs <- probandFilteredCNVs
    # 
    CH_Data[[proband]]$CNVs <- probandCNVs
    CH_Data[[proband]]$Mother.CNV <- motherCNVs
    CH_Data[[proband]]$Father.CNV <- fatherCNVs
    
    lens <- sapply(CH_Data[[proband]], nrow)
    
    ## Get SNVs data
    for(member in c("proband", "Mother", "Father")){
      if(member == "proband"){
        cnv <- "CNVs"
        snvpath <- probandFilepath
        snvout <- "SNVs"
      }else if(member == "Mother"){
        cnv <- "Mother.CNV"
        snvpath <- motherFilepath
        snvout <- "Mother.SNV"
      }else{
        cnv <- "Father.CNV"
        snvpath <- fatherFilepath
        snvout <- "Father.SNV"
      }
      
      snvs <- fread(snvpath)[high_quality == T]
      if(nrow(snvs) > 0 & nrow(CH_Data[[proband]][[cnv]]) > 0){
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
    
    message(sprintf("%s, %s proband SNVs, %s father SNVs, %s mother SNVs", i, 
                    nrow(CH_Data[[proband]][["SNVs"]]), 
                    nrow(CH_Data[[proband]][["Father.SNV"]]), 
                    nrow(CH_Data[[proband]][["Mother.SNV"]])))
  }
  
  return(CH_Data)
}   

print("getting probands data")

# To get data for the comparison of CH event between probands and their transmitting parents
Get_CH_Events(probands = probands_10P,
              CNVdf = CNVs_10Percent,
              SNVfolder = "exonic_splicing") %>%
  write_yaml(.,"SSC_CH_Data_CNV10P_SNV.yaml")

print("getting unaffected siblings data")

# To get data for comparison of CH event between unaffected siblings and their transmitting parents
Get_CH_Events(probands = Unaffected_siblings10P,
              CNVdf = CNVs_10Percent,
              SNVfolder = "exonic_splicing") %>%
  write_yaml(.,"SSC_CH.unaffectedSiblings_Data_CNV10P_SNV.yaml")

print("Done")

