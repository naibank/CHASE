#Required Packages
#BiocManager::install("GenomicRanges")
#BiocManager::install("Repitools")

library(GenomicRanges)
library(Repitools)
library(yaml)
library(data.table)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

Determine_CNV_MAF <- function(CNVs) {
  CNV_MAF <- pmax(CNVs$CGparentalPercFreq_50percRecipOverlap,
                  pmax(CNVs$cnvnIlmXParentalPercFreq_50percRecipOverlap,
                       pmax(CNVs$erdsIlmXParentalPercFreq_50percRecipOverlap,
                            pmax(CNVs$cnvnIlm2ParentalPercFreq_50percRecipOverlap,
                                 pmax(CNVs$erdsIlm2ParentalPercFreq_50percRecipOverlap,
                                      pmax(CNVs$otgCnvnPercFreq_50percRecipOverlap,
                                           pmax(CNVs$otgErdsPercFreq_50percRecipOverlap,
                                                pmax(CNVs$sscErdsPercFreq_50percRecipOverlap,
                                                     CNVs$sscCnvnPercFreq_50percRecipOverlap, na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T)/100
  return (CNV_MAF)
}

Filter_ASD_Associated_CNVs_and_LOFs <- function(df, ds='MSSNG') {
  ls1 <- fread("MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("MSSNG+SSC.CNVs.tsv")
  ls_exclude <- unique(rbind(ls1, ls2))
  
  if (ds == "MSSNG") {
    file <- "MSSNG_metadata.tsv"
  }
  else if (ds == "SSC") {
    file <- "SSC_metadata.tsv"
  }
  
  full_metadata <- fread(file, data.table = F)
  
  familyID_exclude <- full_metadata$`Family ID`[which(full_metadata$`Sample ID` %in% ls_exclude$Sample)]
  sampleID_exclude <- full_metadata$`Sample ID`[which(full_metadata$`Family ID` %in% familyID_exclude)]
  
  df <- df[which(!df$`Sample.ID` %in% sampleID_exclude),]
  
  return (df)
}

Get_CNVs_and_SNVs_Counts <- function(data, ds = "MSSNG", CNV_freq_filter=1, SNV_freq_filter=1,
                                       alt_base_filter_threshold = 0.9, child='proband') {
  
  childSNVs_count <- data.frame()
  childCNVs_count <- data.frame()
  motherSNVs_count <- data.frame()
  motherCNVs_count <- data.frame()
  fatherSNVs_count <- data.frame()
  fatherCNVs_count <- data.frame()
  
  for(i in 1:length(data)){
    childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
    childSNVs <- data.frame(data[[i]]$SNVs)
    
    fatherCNVs <- data.frame(data[[i]]$Father.CNV)
    fatherSNVs <- data.frame(data[[i]]$Father.SNV)
    
    motherCNVs <- data.frame(data[[i]]$Mother.CNV)
    motherSNVs <- data.frame(data[[i]]$Mother.SNV)
    
    #Determine CNV minor allele frequency
    if (length(childCNVs) > 0) {
      childCNVs$cnvFreq <- Determine_CNV_MAF(childCNVs)
      childCNVs <-  childCNVs[which(childCNVs$cnvFreq < CNV_freq_filter), ]
    }
    if (length(fatherCNVs) > 0) {
      fatherCNVs$cnvFreq <- Determine_CNV_MAF(fatherCNVs)
      fatherCNVs <-  fatherCNVs[which(fatherCNVs$cnvFreq < CNV_freq_filter), ]
    }
    if (length(motherCNVs) > 0) {
      motherCNVs$cnvFreq <- Determine_CNV_MAF(motherCNVs)
      motherCNVs <-  motherCNVs[which(motherCNVs$cnvFreq < CNV_freq_filter), ]
    }
    
    childSNVs <- childSNVs[which(childSNVs$freq_max < SNV_freq_filter), ]
    fatherSNVs <- fatherSNVs[which(fatherSNVs$freq_max < SNV_freq_filter), ]
    motherSNVs <- motherSNVs[which(motherSNVs$freq_max < SNV_freq_filter), ]
    
    if (nrow(childCNVs) > 0) {
      childCNVs_c <- nrow(childCNVs)
      child_df <- data.frame(Sample.ID=childCNVs$Sample.ID[1], Count=childCNVs_c)
      childCNVs_count <- rbind(childCNVs_count, child_df)
      
      if(nrow(childSNVs) > 0) {
        childCNVs_g <- GRanges(childCNVs$chrAnn,
                               IRanges(childCNVs$STARTAnn,
                                       childCNVs$ENDAnn), "*")
        
        names(childSNVs)[11] <- "SampleData"
        
        childSNVs.g <- GRanges(childSNVs$CHROM, IRanges(childSNVs$start, childSNVs$end), "*")
        olap <- findOverlaps(childSNVs.g, childCNVs_g)
        childSNVs <- childSNVs[unique(olap@from), ]
        childSNVs <- subset(childSNVs, childSNVs$alt_fraction >= alt_base_filter_threshold)
        childSNVs_c <- nrow(childSNVs)
        childSNV_df <- data.frame(Sample.ID=childSNVs$X.Sample[1], Count=childSNVs_c)
        childSNVs_count <- rbind(childSNVs_count, childSNV_df) 
      }
    }
    
    if (nrow(motherCNVs) > 0) {
      motherCNVs_c <- nrow(motherCNVs)
      mother_df <- data.frame(Sample.ID=motherCNVs$Sample.ID[1], Count=motherCNVs_c)
      motherCNVs_count <- rbind(motherCNVs_count, mother_df)
      if(nrow(motherSNVs) > 0){
        motherCNVs_g <- GRanges(motherCNVs$chrAnn,
                                IRanges(motherCNVs$STARTAnn,
                                        motherCNVs$ENDAnn), "*")
        
        names(motherSNVs)[11] <- "SampleData"
        
        motherSNVs.g <- GRanges(motherSNVs$CHROM, IRanges(motherSNVs$start, motherSNVs$end), "*")
        olap <- findOverlaps(motherSNVs.g, motherCNVs_g)
        motherSNVs <- motherSNVs[unique(olap@from), ]
        motherSNVs <- subset(motherSNVs, motherSNVs$alt_fraction >= alt_base_filter_threshold)
        motherSNVs_c <- nrow(motherSNVs)
        motherSNV_df <- data.frame(Sample.ID=motherSNVs$X.Sample[1], Count=motherSNVs_c)
        motherSNVs_count <- rbind(motherSNVs_count, motherSNV_df) 
      }
    }

    if (nrow(fatherCNVs) > 0) {
      fatherCNVs_c <- nrow(fatherCNVs)
      father_df <- data.frame(Sample.ID=fatherCNVs$Sample.ID[1], Count=fatherCNVs_c)
      fatherCNVs_count <- rbind(fatherCNVs_count, father_df)
      if(nrow(fatherSNVs) > 0){
        fatherCNVs_g <- GRanges(fatherCNVs$chrAnn,
                                IRanges(fatherCNVs$STARTAnn,
                                        fatherCNVs$ENDAnn), "*")
        
        names(fatherSNVs)[11] <- "SampleData"
        
        fatherSNVs.g <- GRanges(fatherSNVs$CHROM, IRanges(fatherSNVs$start, fatherSNVs$end), "*")
        olap <- findOverlaps(fatherSNVs.g, fatherCNVs_g)
        fatherSNVs <- fatherSNVs[unique(olap@from), ]
        fatherSNVs <- subset(fatherSNVs, fatherSNVs$alt_fraction >= alt_base_filter_threshold)
        fatherSNVs_c <- nrow(fatherSNVs)
        fatherSNV_df <- data.frame(Sample.ID=fatherSNVs$X.Sample[1], Count=fatherSNVs_c)
        fatherSNVs_count <- rbind(fatherSNVs_count, fatherSNV_df) 
      }
    }
    message(i)
  }
  
  # Add samples with 0 CNVs from metadata
  if (ds == "MSSNG") {
    file <- "MSSNG_metadata.tsv"
  }
  else if (ds == "SSC") {
    file <- "SSC_metadata.tsv"
  }
  metadata <- fread(file, data.table = F)
  childCNVs_count <- rbind(childCNVs_count, data.frame(Sample.ID=metadata$'Sample ID'[which((!metadata$'Sample ID' %in% childCNVs_count$Sample.ID)
                                                           & metadata$Relation == child)], Count=0))
  motherCNVs_count <- rbind(motherCNVs_count, data.frame(Sample.ID=metadata$'Sample ID'[which((!metadata$'Sample ID'%in% motherCNVs_count$Sample.ID)
                                                                       & metadata$Relation == 'mother')], Count=0))
  fatherCNVs_count <- rbind(fatherCNVs_count, data.frame(Sample.ID=metadata$'Sample ID'[which((!metadata$'Sample ID' %in% fatherCNVs_count$Sample.ID)
                                                                       & metadata$Relation == 'father')], Count=0))
  
  # Filter ASD associated CNVs and LOFs
  childSNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(childSNVs_count, ds)
  childCNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(childCNVs_count, ds)
  motherSNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(motherSNVs_count, ds)
  motherCNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(motherCNVs_count, ds)
  fatherSNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(fatherSNVs_count, ds)
  fatherCNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(fatherCNVs_count, ds)
  
  return (list(childSNVs_count, 
               childCNVs_count, 
               motherSNVs_count, 
               motherCNVs_count,
               fatherSNVs_count,
               fatherCNVs_count))
}

Get_CNVs_Counts <- function(data_list, ds = "MSSNG", CNV_freq_filter=0.1, SNV_freq_filter=1,
                                    alt_base_filter_threshold = 0.9, child='proband') {
  
  childCNVs_count <- data.frame()
  motherCNVs_count <- data.frame()
  fatherCNVs_count <- data.frame()
  for (j in 1:length(data_list)) {
    data <- data_list[[j]]
    for(i in 1:length(data)){
      childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
      
      fatherCNVs <- data.frame(data[[i]]$Father.CNV)
      
      motherCNVs <- data.frame(data[[i]]$Mother.CNV)
      
      #Determine CNV minor allele frequency
      if (length(childCNVs) > 0) {
        childCNVs$cnvFreq <- Determine_CNV_MAF(childCNVs)
        childCNVs <-  childCNVs[which(childCNVs$cnvFreq <= CNV_freq_filter), ]
      }
      if (length(fatherCNVs) > 0) {
        fatherCNVs$cnvFreq <- Determine_CNV_MAF(fatherCNVs)
        fatherCNVs <-  fatherCNVs[which(fatherCNVs$cnvFreq <= CNV_freq_filter), ]
      }
      if (length(motherCNVs) > 0) {
        motherCNVs$cnvFreq <- Determine_CNV_MAF(motherCNVs)
        motherCNVs <-  motherCNVs[which(motherCNVs$cnvFreq <= CNV_freq_filter), ]
      }
      
      if (nrow(childCNVs) > 0) {
        childCNVs_c <- nrow(childCNVs)
        child_df <- data.frame(Sample.ID=childCNVs$Sample.ID[1], Count=childCNVs_c)
        childCNVs_count <- rbind(childCNVs_count, child_df)
      }
      
      if (nrow(motherCNVs) > 0) {
        motherCNVs_c <- nrow(motherCNVs)
        mother_df <- data.frame(Sample.ID=motherCNVs$Sample.ID[1], Count=motherCNVs_c)
        motherCNVs_count <- rbind(motherCNVs_count, mother_df)
      }
      
      if (nrow(fatherCNVs) > 0) {
        fatherCNVs_c <- nrow(fatherCNVs)
        father_df <- data.frame(Sample.ID=fatherCNVs$Sample.ID[1], Count=fatherCNVs_c)
        fatherCNVs_count <- rbind(fatherCNVs_count, father_df)
      }
      message(i)
    }
  }
  # Filter ASD associated CNVs and LOFs
  childCNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(childCNVs_count, ds)
  motherCNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(motherCNVs_count, ds)
  fatherCNVs_count <- Filter_ASD_Associated_CNVs_and_LOFs(fatherCNVs_count, ds)
  
  return (list(childCNVs_count, 
               motherCNVs_count,
               fatherCNVs_count))
}

Save_Boxplot <- function(filename, form, df, main, notch=T, names, outline=T) {
  png(filename)
  boxplot(form, data=df, main=main, notch=notch, names=names, outline=outline)
  dev.off()
}

#### Get raw CNV counts ####
ILMN_data <- yaml::yaml.load_file("./Unbiased/MSSNG_ILMN_CH_Data_CNV10P_SNV.yaml")
CG_data <- yaml::yaml.load_file("./Unbiased/MSSNG_CG_CH_Data_CNV10P_SNV.yaml")
MSSNG_data_list <- list(ILMN_data, CG_data)
MSSNG_CNV_Counts <- Get_CNVs_Counts(data_list=MSSNG_data_list)
MSSNG_proband_counts <- MSSNG_CNV_Counts[[1]]
MSSNG_mother_counts <- MSSNG_CNV_Counts[[2]]
MSSNG_father_counts <- MSSNG_CNV_Counts[[3]]

SSC_data_list <- yaml::yaml.load_file("./Unbiased/SSC_CH_Data_CNV10P_SNV.yaml")
SSC_CNV_Counts <- Get_CNVs_Counts(data_list=list(SSC_data_list), ds='SSC')
SSC_proband_counts <- SSC_CNV_Counts[[1]]
SSC_mother_counts <- SSC_CNV_Counts[[2]]
SSC_father_counts <- SSC_CNV_Counts[[3]]

SSC_US_data_list <- yaml::yaml.load_file("./Unbiased/SSC_CH.unaffectedSiblings_Data_CNV10P_SNV.yaml")
SSC_US_CNV_Counts <- Get_CNVs_Counts(data_list=list(SSC_US_data_list), ds='SSC')
SSC_US_counts <- SSC_US_CNV_Counts[[1]]
SSC_US_mother_counts <- SSC_US_CNV_Counts[[2]]
SSC_US_father_counts <- SSC_US_CNV_Counts[[3]]

##### ILMN DATA (1111 CH events) ##### 
ILMN_data <- yaml::yaml.load_file("./Unbiased/MSSNG_ILMN_CH_Data_CNV10P_SNV.yaml")
ILMN_CNV_SNV_data <- Get_CNVs_and_SNVs_Counts(ILMN_data, ds="MSSNG", CNV_freq_filter=0.01, SNV_freq_filter=1)
ILMN_proband_SNVs_count <- ILMN_CNV_SNV_data[[1]]
ILMN_proband_CNVs_count <- ILMN_CNV_SNV_data[[2]]
ILMN_mother_SNVs_count <- ILMN_CNV_SNV_data[[3]]
ILMN_mother_CNVs_count <- ILMN_CNV_SNV_data[[4]]
ILMN_father_SNVs_count <- ILMN_CNV_SNV_data[[5]]
ILMN_father_CNVs_count <- ILMN_CNV_SNV_data[[6]]

# Convert to gender and stack
# ILMN
colnames(ILMN_proband_CNVs_count)[which(names(ILMN_proband_CNVs_count) == 'Sample.ID')] <- 'Relation'
ILMN_proband_CNVs_count$Relation <- 0
colnames(ILMN_mother_CNVs_count)[which(names(ILMN_mother_CNVs_count) == 'Sample.ID')] <- 'Relation'
ILMN_mother_CNVs_count$Relation <- 1
colnames(ILMN_father_CNVs_count)[which(names(ILMN_father_CNVs_count) == 'Sample.ID')] <- 'Relation'
ILMN_father_CNVs_count$Relation <- 2
ILMN_CNVs_count_by_relation <- rbind(ILMN_proband_CNVs_count, ILMN_mother_CNVs_count, ILMN_father_CNVs_count)

colnames(ILMN_proband_SNVs_count)[which(names(ILMN_proband_SNVs_count) == 'Sample.ID')] <- 'Relation'
ILMN_proband_SNVs_count$Relation <- 0
colnames(ILMN_mother_SNVs_count)[which(names(ILMN_mother_SNVs_count) == 'Sample.ID')] <- 'Relation'
ILMN_mother_SNVs_count$Relation <- 1
colnames(ILMN_father_SNVs_count)[which(names(ILMN_father_SNVs_count) == 'Sample.ID')] <- 'Relation'
ILMN_father_SNVs_count$Relation <- 2
ILMN_SNVs_count_by_relation <- rbind(ILMN_proband_SNVs_count, ILMN_mother_SNVs_count, ILMN_father_SNVs_count)

Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/ILMN_CNVs_count_by_relation.png", Count~Relation, ILMN_CNVs_count_by_relation,
             main='ILMN CNV Counts by Relation', names = c('proband', 'mother', 'father'))
Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/ILMN_SNVs_count_by_relation.png", Count~Relation, ILMN_SNVs_count_by_relation,
             main='ILMN SNV Counts by Relation', names = c('proband', 'mother', 'father'))

##### CG DATA ##### (1111 CH events)
CG_data <- yaml::yaml.load_file("./Unbiased/MSSNG_CG_CH_Data_CNV10P_SNV.yaml")
CG_CNV_SNV_data <- Get_CNVs_and_SNVs_Counts(CG_data, ds="MSSNG", CNV_freq_filter=0.01, SNV_freq_filter=1)
CG_proband_SNVs_count <- CG_CNV_SNV_data[[1]]
CG_proband_CNVs_count <- CG_CNV_SNV_data[[2]]
CG_mother_SNVs_count <- CG_CNV_SNV_data[[3]]
CG_mother_CNVs_count <- CG_CNV_SNV_data[[4]]
CG_father_SNVs_count <- CG_CNV_SNV_data[[5]]
CG_father_CNVs_count <- CG_CNV_SNV_data[[6]]

# Convert to gender and stack
# CG
colnames(CG_proband_CNVs_count)[which(names(CG_proband_CNVs_count) == 'Sample.ID')] <- 'Relation'
CG_proband_CNVs_count$Relation <- 0
colnames(CG_mother_CNVs_count)[which(names(CG_mother_CNVs_count) == 'Sample.ID')] <- 'Relation'
CG_mother_CNVs_count$Relation <- 1
colnames(CG_father_CNVs_count)[which(names(CG_father_CNVs_count) == 'Sample.ID')] <- 'Relation'
CG_father_CNVs_count$Relation <- 2
CG_CNVs_count_by_relation <- rbind(CG_proband_CNVs_count, CG_mother_CNVs_count, CG_father_CNVs_count)

colnames(CG_proband_SNVs_count)[which(names(CG_proband_SNVs_count) == 'Sample.ID')] <- 'Relation'
CG_proband_SNVs_count$Relation <- 0
colnames(CG_mother_SNVs_count)[which(names(CG_mother_SNVs_count) == 'Sample.ID')] <- 'Relation'
CG_mother_SNVs_count$Relation <- 1
colnames(CG_father_SNVs_count)[which(names(CG_father_SNVs_count) == 'Sample.ID')] <- 'Relation'
CG_father_SNVs_count$Relation <- 2
CG_SNVs_count_by_relation <- rbind(CG_proband_SNVs_count, CG_mother_SNVs_count, CG_father_SNVs_count)

Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/CG_CNVs_count_by_relation.png", Count~Relation, CG_CNVs_count_by_relation,
             main='CG CNV Counts by Relation', names = c('proband', 'mother', 'father'))
Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/CG_SNVs_count_by_relation.png", Count~Relation, CG_SNVs_count_by_relation,
             main='CG SNV Counts by Relation', names = c('proband', 'mother', 'father'))

##### SSC DATA ##### (1111 CH events)
SSC_data <- yaml::yaml.load_file("./Unbiased/SSC_CH_Data_CNV10P_SNV.yaml")
SSC_CNV_SNV_data <- Get_CNVs_and_SNVs_Counts(SSC_data, ds="SSC", CNV_freq_filter=0.01, SNV_freq_filter=1)
SSC_proband_SNVs_count <- SSC_CNV_SNV_data[[1]]
SSC_proband_CNVs_count <- SSC_CNV_SNV_data[[2]]
SSC_mother_SNVs_count <- SSC_CNV_SNV_data[[3]]
SSC_mother_CNVs_count <- SSC_CNV_SNV_data[[4]]
SSC_father_SNVs_count <- SSC_CNV_SNV_data[[5]]
SSC_father_CNVs_count <- SSC_CNV_SNV_data[[6]]

# Convert to gender and stack
# SSC
colnames(SSC_proband_CNVs_count)[which(names(SSC_proband_CNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_proband_CNVs_count$Relation <- 0
colnames(SSC_mother_CNVs_count)[which(names(SSC_mother_CNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_mother_CNVs_count$Relation <- 1
colnames(SSC_father_CNVs_count)[which(names(SSC_father_CNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_father_CNVs_count$Relation <- 2
SSC_CNVs_count_by_relation <- rbind(SSC_proband_CNVs_count, SSC_mother_CNVs_count, SSC_father_CNVs_count)

colnames(SSC_proband_SNVs_count)[which(names(SSC_proband_SNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_proband_SNVs_count$Relation <- 0
colnames(SSC_mother_SNVs_count)[which(names(SSC_mother_SNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_mother_SNVs_count$Relation <- 1
colnames(SSC_father_SNVs_count)[which(names(SSC_father_SNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_father_SNVs_count$Relation <- 2
SSC_SNVs_count_by_relation <- rbind(SSC_proband_SNVs_count, SSC_mother_SNVs_count, SSC_father_SNVs_count)

Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/SSC_CNVs_count_by_relation.png", Count~Relation, SSC_CNVs_count_by_relation,
             main='SSC CNV Counts by Relation', names = c('proband', 'mother', 'father'))
Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/SSC_SNVs_count_by_relation.png", Count~Relation, SSC_SNVs_count_by_relation,
             main='SSC SNV Counts by Relation', names = c('proband', 'mother', 'father'))

##### SSC unaffected sibs DATA ##### (1111 CH events)
SSC_data <- yaml::yaml.load_file("./Unbiased/SSC_CH.unaffectedSibilngs_Data_CNV10P_SNV.yaml")
SSC_CNV_SNV_data <- Get_CNVs_and_SNVs_Counts(SSC_data, ds="SSC", CNV_freq_filter=0.01, SNV_freq_filter=1,
                                             child='unaffected sibling')
SSC_US_SNVs_count <- SSC_CNV_SNV_data[[1]]
SSC_US_CNVs_count <- SSC_CNV_SNV_data[[2]]
SSC_motherUS_SNVs_count <- SSC_CNV_SNV_data[[3]]
SSC_motherUS_CNVs_count <- SSC_CNV_SNV_data[[4]]
SSC_fatherUS_SNVs_count <- SSC_CNV_SNV_data[[5]]
SSC_fatherUS_CNVs_count <- SSC_CNV_SNV_data[[6]]

# Convert to gender and stack
# SSC
colnames(SSC_US_CNVs_count)[which(names(SSC_US_CNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_US_CNVs_count$Relation <- 0
colnames(SSC_motherUS_CNVs_count)[which(names(SSC_motherUS_CNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_motherUS_CNVs_count$Relation <- 1
colnames(SSC_fatherUS_CNVs_count)[which(names(SSC_fatherUS_CNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_fatherUS_CNVs_count$Relation <- 2
SSC_USCNVs_count_by_relation <- rbind(SSC_US_CNVs_count, SSC_motherUS_CNVs_count, SSC_fatherUS_CNVs_count)

colnames(SSC_US_SNVs_count)[which(names(SSC_US_SNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_US_SNVs_count$Relation <- 0
colnames(SSC_motherUS_SNVs_count)[which(names(SSC_motherUS_SNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_motherUS_SNVs_count$Relation <- 1
colnames(SSC_fatherUS_SNVs_count)[which(names(SSC_fatherUS_SNVs_count) == 'Sample.ID')] <- 'Relation'
SSC_fatherUS_SNVs_count$Relation <- 2
SSC_USSNVs_count_by_relation <- rbind(SSC_US_SNVs_count, SSC_motherUS_SNVs_count, SSC_fatherUS_SNVs_count)

Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/SSC_USCNVs_count_by_relation.png", Count~Relation, SSC_USCNVs_count_by_relation,
             main='SSC (Unaffected Sibs) CNV Counts by Relation', names = c('unaffected sibling', 'mother', 'father'))
Save_Boxplot("../DT/Unbiased/Plots/zero_cnv/SSC_USSNVs_count_by_relation.png", Count~Relation, SSC_USSNVs_count_by_relation,
             main='SSC (Unaffected Sibs) SNV Counts by Relation', names = c('unaffected sibling', 'mother', 'father'))
