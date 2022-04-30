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

# setwd("/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG")

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

# Getting CNVs (with different MAF limits) for samples sequenced by ILMN (only CNVs_10Percent is used further on)
# CNVs_01Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG/CNVs.ILMN/CNVs.MSSNG_ILMN.freq_0.1percent.HQR.tsv")
# CNVs_1Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG/CNVs.ILMN/CNVs.MSSNG_ILMN.freq_1percent.HQR.tsv")
# CNVs_5Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG/CNVs.ILMN/CNVs.MSSNG_ILMN.freq_5percent.HQR.tsv")
CNVs_10Percent <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/CNVs.ILMN/CNVs.MSSNG_ILMN.freq_10percent.HQR.IntFreq.tsv")

# Getting CNVs (with different MAF limits) for samples sequenced by CG platform
CNVs_10Percent.CG <- CNVfilter.procedure(CNVfile.path = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/CNVs.CG/CNVs.MSSNG_CG.freq_10percent.HQR.IntFreq.tsv")

## STEP 2: proband filtering 
# Exclude families when proband also has an ASD affected parent
Family_to_exclude <- as_tibble(read.delim("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/MSSNG_metadata.tsv", stringsAsFactors = F)) %>%
  filter(Relation == "proband") %>%
  filter(Sample.ID %in% Mother.ID | Sample.ID %in% Father.ID) %>% pull(Family.ID)

# Subjects to exclude that failed CNV quality control (list that Brett supplied)
Subjects_to_exclude <- c("1-0004-004","1-0059-001A","1-0075-003","1-0079-004"," 1-0092-003","1-0202-003"," 1-0219-004",
                         "1-0222-003","1-0229-001","1-0292-003","1-0298-002","1-0337-002","1-0360-003"," 1-0366-003","1-0385-001",
                         "1-0389-003","1-0394-002","1-0401-003","1-0414-003","1-0414-004","1-0414-006","1-0432-003","1-0494-004","1-0534-001","1-0534-003","1-0534-006","1-0562-004","1-0574-001",
                         "1-0618-002","1-0629-002","1-0638-003","1-0647-001","1-0658-003","1-0680-001","1-0695-001","1-0714-002","1-0731-004",
                         "1-0816-004","1-0838-001","1-0844-001","1-0844-002","1-0919-003A","1-0960-002","1-0974-003","1-0998-002","1-1015-001","1-1016-002","10-0011-003",
                         "10-0013-001","10-1120-002","10-1127-001","2-0019-003","2-0142-003","2-0143-005","2-0160-002","2-0198-002",
                         "2-0215-004","2-0241-003","2-0319-001","2-1117-001","2-1139-002","2-1139-003","2-1148-003","2-1172-001","2-1205-001","2-1210-003",
                         "2-1220-002","2-1220-003","2-1245-001","2-1245-003","2-1265-001","2-1272-001A","2-1272-002A","2-1275-003",
                         "2-1302-001","2-1313-003","2-1315-001","2-1335-003","2-1355-001","2-1381-002","2-1438-003","2-1475-001","2-1502-003","2-1514-002","2-1592-002",
                         "2-1592-003","2-1723-003","2-1727-004","3-0135-000","3-0141-000","3-0141-100","3-0168-000","3-0204-000","3-0216-000","3-0612-000","3-0666-102",
                         "3-0791-100","4-0021-004","4-0042-001","4-0062-003","4-0070-004","4-0073-003","4-0073-005","4-0079-001","4-0082-004","5-0005-001","5-0009-003",
                         "5-0111-001","5-5004-002","5-5072-003","7-0068-003","7-0097-003","7-0104-002","7-0123-001","7-0151-003","7-0232-004","7-0277-001","7-0277-002",
                         "7-0277-003","7-0286-001","7-0295-002","7-0314-001","7-0322-003","AM01ZF-02","AU0039201","AU0039202","AU004403","AU023012","AU055603","AU059003",
                         "AU073001","AU1635301","AU2137201","AU2142301","AU2218201","AU2218202","AU2218301","AU2272301","AU2283301","AU2295301","AU2300301","AU2302301",
                         "AU2303301","AU2310301","AU2320301","AU2385301","AU2410302","AU2463301","AU3610201","AU3695201","AU3809202","AU4145301","MSSNG00012-003","MSSNG00019-002",
                         "MSSNG00024-001","MSSNG00024-002","MSSNG00024-003A","MSSNG00047-002","MSSNG00384-002","MSSNG00395-003","MSSNG00403-003","MSSNG00432-003","MT_107.3",
                         "MT_12.2","MT_121.1","MT_121.3","MT_13.3","MT_165.2","MT_172.1","MT_172.2","MT_182.2","MT_67.3","REACH000315","REACH000322","REACH000336","REACH000337",
                         "REACH000340","REACH000341","REACH000342","REACH000343","REACH000344","REACH000345","REACH000346","REACH000347","REACH000348","REACH000349","REACH000350",
                         "REACH000352","REACH000354","REACH000355","REACH000356","REACH000357","REACH000366","REACH000367","REACH000399","REACH000435","REACH000436","REACH000437",
                         "REACH000442","REACH000443","REACH000448","REACH000449","REACH000455","REACH000457","REACH000467","REACH000468","REACH000471","REACH000475","REACH000476",
                         "REACH000477","REACH000478","REACH000479","REACH000480","REACH000483","REACH000484","REACH000486","REACH000487","REACH000488","REACH000490","REACH000496",
                         "REACH000501","REACH000502","REACH000503","REACH000504","REACH000505","REACH000507","REACH000510","REACH000511","REACH000513","REACH000514","REACH000515","REACH000516",
                         "REACH000517","REACH000518","REACH000519","REACH000520","REACH000521","REACH000523","REACH000524","REACH000525","REACH000526","REACH000528","REACH000529","REACH000530","REACH000555",
                         "REACH000556","REACH000557","REACH000558","REACH000559","REACH000560","REACH000561","REACH000562","REACH000563","REACH000565","REACH000578","REACH000580","REACH000584","REACH000588","REACH000592","REACH000593","REACH000597","REACH000598",
                         "REACH000601","REACH000604","REACH000610","REACH000612","REACH000626","REACH000631","REACH000632","REACH000633","REACH000639","REACH000648","REACH000651",
                         "REACH000652","REACH000654","REACH000656","REACH000662","REACH000666","REACH000671","REACH000672","REACH000674","REACH000675","REACH000676","REACH000677",
                         "REACH000679","REACH000680","REACH000683","REACH000688","REACH000693","REACH000701","REACH000702","REACH000705","REACH000706","REACH000708","REACH000710","REACH000714",
                         "REACH000715","REACH000716","REACH000717","REACH000721","REACH000722","REACH000760","REACH000764","REACH000768","REACH000769","REACH000770","REACH000772",
                         "REACH000773","REACH000774","REACH000775","REACH000776","SJD_34.4","T2T9E-03",'1-0142-001','1-0298-002','1-0321-004','1-0401-003',
                         '1-0534-001','1-0534-003','1-0534-006','2-0142-003','2-0143-005','2-0198-002','2-0319-001','2-1148-003','2-1210-003','2-1252-003','2-1275-003','2-1335-003','2-1355-001','2-1723-003')

# Function to locate probands that have at least one CNV that matches the criteria in the CNV filter procedure function
# Insertion of 3 parameters are necessary:
# Metadata.filepath = path to metadata file
# CNV.database = database of CNVs that remain after filter procedure. The rest is self explanatory 
LocateProbands <- function(metadata.filepath, CNV.database, Family_to_exclude, Subjects_to_exclude){
  # Selects all probands, filters out samples that are re-sequenced and that have to be excluded 
  probands <- as_tibble(read.delim(metadata.filepath, stringsAsFactors = FALSE)) %>%
    dplyr::filter(Relation == "proband") %>%
    filter(Exclude.because.re.sequenced. != "yes") %>%
    filter(!Family.ID %in% Family_to_exclude &
             !Sample.ID %in% Subjects_to_exclude) #Exclude probands if their sample.ID is located in the vector of individuals that need to be excluded
  
  # Only keeps probands who have at least one CNV in the CNV database after filtering procedure 
  probands <- probands %>%
    dplyr::filter(Sample.ID %in% CNV.database$Sample.ID)
  return(probands)
}

print("TEST 2")

## Get probands using 10 Percent CNV database sequenced by ILMN
probands_10P <- LocateProbands(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/MSSNG_metadata.tsv",
                                 CNV.database = CNVs_10Percent,
                                 Family_to_exclude = Family_to_exclude,
                                 Subjects_to_exclude = Subjects_to_exclude)
# write_xlsx(probands_10P, "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG//probands_10P.xlsx")

## Get probands using 10 percent CNV database sequenced by CG
probands_10P.CG <- LocateProbands(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/MSSNG_metadata.tsv",
                                 CNV.database = CNVs_10Percent.CG,
                                 Family_to_exclude = Family_to_exclude,
                                 Subjects_to_exclude = Subjects_to_exclude)
# write_xlsx(probands_10P.CG, "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG//probands_10P_CG.xlsx")

## Locate parent (IDs) of probands that have at least one inherited CNV through Family.ID column
LocateParents <- function(metadata.filepath, probands.data, Family_to_exclude, Subjects_to_exclude){
  ## Filters out parents/samples that have been resequenced, are not in the vectors of samples that need to be excluded and keeps only the parents of the probands with at least one inherited CNV deletion
  parents <- as_tibble(read.delim(metadata.filepath, stringsAsFactors = FALSE)) %>%
  filter(Exclude.because.re.sequenced. != "yes") %>%
    filter(!Family.ID %in% Family_to_exclude &
             !Sample.ID %in% Subjects_to_exclude) %>%
    dplyr::filter(Sample.ID %in% probands.data$Father.ID | Sample.ID %in% probands.data$Mother.ID)
  
  return(parents)
}

# Get samples of probands' parents sequenced by ILMN platform
parents_10P <- LocateParents(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/MSSNG_metadata.tsv",
                               probands.data = probands_10P,
                               Family_to_exclude = Family_to_exclude,
                               Subjects_to_exclude = Subjects_to_exclude)
#write_xlsx(parents_10P, "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG//parents_10P.xlsx")

# Get samples of probands' parents sequenced by CG platform
parents_10P.CG <- LocateParents(metadata.filepath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/MSSNG_metadata.tsv",
                               probands.data = probands_10P.CG,
                               Family_to_exclude = Family_to_exclude,
                               Subjects_to_exclude = Subjects_to_exclude)
#write_xlsx(parents_10P.CG, "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG//parents_10P_CG.xlsx")

print("TEST 3")

#################################################################################################################################################

## STEP 3a: Rerun analysis with pairs (mother + proband, father + proband) instead of trios
# Need to filter out probands that don't have both parents sequenced
# Pull family ID from all parents of probands that have at least one inherited CNV (parents_10p and parents_10P.CG)
bothParents <- function(parents){
  parentID <- vector("character", length = nrow(parents))
  for(i in 1:nrow(parents)){
    parentID[i] <- parents$Family.ID[i]
  }
  return(parentID)
}

bothParents_10P <- bothParents(parents = parents_10P)
bothParents_10P.CG <- bothParents(parents = parents_10P.CG)

# Sort family ID of parents to match identical IDs (matched IDs indicate both parents are sequenced)
bothParents_10P <- sort.int(bothParents_10P, partial = NULL, decreasing = FALSE, method = c("auto"), index.return = FALSE)
bothParents_10P.CG <- sort.int(bothParents_10P.CG, partial = NULL, decreasing = FALSE, method = c("auto"), index.return = FALSE)

# Returns logical matrix ("FALSE" = no duplicate, "TRUE" = duplicate i.e. both parents sequenced )
parentDuplicates_10P <- duplicated(bothParents_10P, incomparables = FALSE)
parentDuplicates_10P.CG <- duplicated(bothParents_10P.CG, incomparables = FALSE)

# Filter bothParents to obtain Family IDs of probands with both parents sequenced
filterParents <- function(duplicates, bothParents){
  for(i in 1:NROW(bothParents)){ # NROW for vectors
    if(duplicates[i] == "FALSE"){
      bothParents[i] = NA # Replaces Family IDs of single parents with "NA" in bothParents
    }
  }
  bothParents = bothParents[ !is.na(bothParents) ] #Drop NAs (left with Family IDs of parent pairs)
  return(bothParents)
}

filteredParents_10P <- data.frame(filterParents(duplicates = parentDuplicates_10P, bothParents = bothParents_10P))
filteredParents_10P.CG <- data.frame(filterParents(duplicates = parentDuplicates_10P.CG, bothParents = bothParents_10P.CG))

names(filteredParents_10P)[1] <- 'Family.ID' #Change column name to match probands to use intersect() function later
names(filteredParents_10P.CG)[1] <- 'Family.ID'

# Pair probands from 10% CNV database with parent doubles (i.e. mother + father) using Family ID
familyTrio <- function(probands, filteredParents){
  
  completeFamily <- vector("character", length = 1)
  completeFamily <- intersect(probands$Family.ID, filteredParents$Family.ID)
  
  return(completeFamily)
}

familyTrio_10P <- data.frame(familyTrio(probands = probands_10P, filteredParents = filteredParents_10P))
familyTrio_10P.CG <- data.frame(familyTrio(probands = probands_10P.CG, filteredParents = filteredParents_10P.CG))

names(familyTrio_10P)[1] <- 'Family.ID'
names(familyTrio_10P.CG)[1] <- 'Family.ID'

# Filter out probands who do not have both parents sequenced
probandsRemoved <- function(probands, familyTrio){
  removeProbands <- probands$Family.ID %in% familyTrio$Family.ID #returns a logical matrix where "TRUE" indicates the probands' family ID is included in the familyTrio
  probands <- subset(probands, removeProbands == "TRUE") #take all 'TRUE' values and store in "probands"
  
  return(probands)
}

probandsRemoved_10P <- probandsRemoved(probands = probands_10P, familyTrio = familyTrio_10P)
probandsRemoved_10P.CG <- probandsRemoved(probands = probands_10P.CG, familyTrio = familyTrio_10P.CG)

# Probands in probandsRemoved data frames all have CNV deletion (and both mom + dad sequenced), now find probands with SNVs to determine occurrence of CH events


Get_CH_Events <- function(probands, CNVdf, SNVfolder){
  CH_Data <- list()
  write.csv(probands$Sample.ID, "TestSampleIDs.csv")
  for(i in 1:nrow(probands)){
    
    # Get the filepath to the proband's relevant SNV file
    probandFilepath <- str_c(c("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/", SNVfolder, "/", 
                               probands$Sample.ID[i], ".tsv.gz"), collapse = "") # Search proband SNV data
                             
    proband <- probands$Sample.ID[i]
    
    # Empty vector for probands' mother
    mother <- vector("character", length = 1)
    
    # Get the filepath to mother SNV file
    mother <- probands %>% dplyr::filter(probands$Sample.ID == probands$Sample.ID[i]) %>% 
      pull(Mother.ID) # Search probands' mother SNV data
    motherFilepath <- str_c(c("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/", SNVfolder, "/",
                              mother[1], ".tsv.gz"), collapse = "")
    
    # Empty vector for probands' father
    father <- vector("character", length = 1)
    
    # Get the filepath to father SNV file
    father <- probands %>% dplyr::filter(probands$Sample.ID == probands$Sample.ID[i]) %>%
      pull(Father.ID) #search probands' father SNV data
    fatherFilepath <- str_c(c("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/", SNVfolder, "/",
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
  
    CH_Data[[proband]]$CNVs <- probandCNVs
    CH_Data[[proband]]$Mother.CNV <- motherCNVs
    CH_Data[[proband]]$Father.CNV <- fatherCNVs
    
    lens <- sapply(CH_Data[[proband]], nrow)
    
    # Get SNVs data
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
      
      tryCatch(
        expr = {
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
        },
        error = function(e){ 
          message(paste(c("Error reading file:", snvpath), collapse=" "))
          message(e)
        },
        warning = function(w){
          message(paste(c("Warning reading file:", snvpath), collapse=" "))
          message(w)
        },
        finally = {
          next
        }
      )
    }
    
    message(sprintf("%s, %s proband SNVs, %s father SNVs, %s mother SNVs", i, 
                    nrow(CH_Data[[proband]][["SNVs"]]), 
                    nrow(CH_Data[[proband]][["Father.SNV"]]), 
                    nrow(CH_Data[[proband]][["Mother.SNV"]])))
  }
  
  return(CH_Data)
}   

# 
# TEST_10P <- getCH_events(probands = probandsRemoved_10P, 
#                          CNVdf = CNVs_10Percent, 
#                          CNVpath = "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG/CNVs.ILMN/CNVs.MSSNG_ILMN.freq_10percent.HQR.tsv", 
#                          SNVfolder = "exonic_splicing.ILMN")

##ILMN

print("getting ILMN probands data")
Get_CH_Events(probands = probands_10P,
             CNVdf = CNVs_10Percent,
             SNVfolder = "exonic_splicing.ILMN") %>%
 write_yaml(.,"MSSNG_ILMN_CH_Data_CNV10P_SNV.yaml")

## CG

print("getting CG probands data")
Get_CH_Events(probands = probands_10P.CG,
              CNVdf = CNVs_10Percent.CG,
              SNVfolder = "exonic_splicing.CG") %>%
  write_yaml(.,"MSSNG_CG_CH_Data_CNV10P_SNV.yaml")

print("Done")

