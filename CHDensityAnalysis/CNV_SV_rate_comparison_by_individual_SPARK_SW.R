################################################################################################################################################################################
# 
# CNV_SV_rate_comparison_by_individual_SW_SPARK.R
# purpose: - creates CH_hits tables (CH hit and total deleted exonic size tables)
#          - plots SPARK linear del size vs. CH count plots
#          - outputs SPARK Fisher's test results using linear observed line as cut-off 
# input: SPARK SPARK_WGS_*1-3*_metadata_relfixed.tsv, CRVs, excluded.qcfailed.samples.txt, SPARKWGS*1-3*_parent_proband_SNVsExonicSizes.yaml
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis
# output: SPARK del size vs. CH count plots & Fisher's test results using linear 
#           del size vs. CH count as cut-off 
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/figures/
# 
# notes: 
#
##############################################################################################################################################################################

library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)

setwd("/Users/shaniawu/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/")

Get_Failed_QC_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_BA/data/excluded.qcfailed.samples.txt", header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) 
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/SPARK/data_processing/data/samples.with.crCNV.SPARK.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) 
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/Data/MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/Data/MSSNG+SSC.CNVs.tsv")
  ls_exclude <- unique(rbind(ls1, ls2))
  
  full_metadata <- fread(file, data.table = F)
  familyID_exclude <- full_metadata$`Family ID`[which(full_metadata$`Sample ID` %in% ls_exclude$Sample)]
  sampleID_exclude <- full_metadata$`Sample ID`[which(full_metadata$`Family ID` %in% familyID_exclude)]
  
  return (sampleID_exclude)
}

Filter_Metadata <- function(file) {
  metadata <- read.delim(file, stringsAsFactors = F)
  
  # Get samples that failed QC
  failed_QC <- Get_Failed_QC_Samples()
  # Get samples with cr CNV
  samples_with_cr_CNV <- Get_cr_CNV_Samples()
  # Get families to exclude that have parents AS the proband as well
  families_to_exclude <- metadata$Family.ID[(metadata$Sample.ID[metadata$Relation == 'proband'] %in% metadata$Mother.ID) |
                                              (metadata$Sample.ID[metadata$Relation == 'proband'] %in% metadata$Father.ID)]
  # Get samples with ASD associated CNVs and LOFs
  # ASD_associated_samples <- Get_SampleID_ASD_Associated_CNVs_and_LOFs(file)
  
  # Filter out from metadata
  metadata <- metadata[!metadata$Sample.ID %in% failed_QC, ]
  metadata <- metadata[!metadata$Sample.ID %in% samples_with_cr_CNV, ]
  metadata <- metadata[!metadata$Family.ID %in% families_to_exclude, ]
  # metadata <- metadata[!metadata$Sample.ID %in% ASD_associated_samples, ]
  
  # List of samples belonging to families with a proband and both parents data:
  probands <- metadata$Family.ID[which(metadata$Relation %in% c("proband", "affected sibling"))]
  US <- metadata$Family.ID[which(metadata$Relation %in% c("unaffected sibling", "other sibling"))]
  mothers <- metadata$Family.ID[which(metadata$Relation == "mother")]
  fathers <- metadata$Family.ID[which(metadata$Relation == "father")]  
  proband_families <- intersect(intersect(probands, mothers), fathers)  # a vector of Family IDs with a full family dataset
  US_families <- intersect(intersect(US, mothers), fathers)  # a vector of Family IDs with a full family dataset
  
  metadata <- metadata[which(metadata$Family.ID %in% c(proband_families,US_families)), ]
  
  return (metadata)
}

### Get df of SNVs and ExonicSizes:
# WGS1
SPARK_parent_proband_SNVsExonicSizes1 <- yaml::yaml.load_file("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/SPARKWGS1_parent_proband_SNVsExonicSizes.yaml")
SPARK_parent_proband_SNVs1 <- data.frame(SPARK_parent_proband_SNVsExonicSizes1[[1]])
SPARK_parent_proband_SNVs1 <- SPARK_parent_proband_SNVs1[SPARK_parent_proband_SNVs1$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')## filerting SNVs?
                                                       | SPARK_parent_proband_SNVs1$LoF,] 
SPARK_parent_proband_SNVs1 <- SPARK_parent_proband_SNVs1[SPARK_parent_proband_SNVs1$cnvENDAnn - SPARK_parent_proband_SNVs1$cnvSTARTAnn <= 10000, ] # only keep SNVs in CNVs <= 10,000 bp
SPARK_parent_proband_SNVs1$UID <- paste(SPARK_parent_proband_SNVs1$X.Sample, SPARK_parent_proband_SNVs1$X.id, sep='.')
SPARK_parent_proband_SNVs1 <- SPARK_parent_proband_SNVs1[!(SPARK_parent_proband_SNVs1$LoF & SPARK_parent_proband_SNVs1$freq_max > 0.01),] #remove LoF with snv freq > 1%

SPARK_parent_proband_ExonicSizes1 <- data.frame(SPARK_parent_proband_SNVsExonicSizes1[[2]])

# WGS2
SPARK_parent_proband_SNVsExonicSizes2 <- yaml::yaml.load_file("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/SPARKWGS2_parent_proband_SNVsExonicSizes.yaml")
SPARK_parent_proband_SNVs2 <- data.frame(SPARK_parent_proband_SNVsExonicSizes2[[1]])
SPARK_parent_proband_SNVs2 <- SPARK_parent_proband_SNVs2[SPARK_parent_proband_SNVs2$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')## filerting SNVs?
                                                         | SPARK_parent_proband_SNVs2$LoF,] 
SPARK_parent_proband_SNVs2 <- SPARK_parent_proband_SNVs2[SPARK_parent_proband_SNVs2$cnvENDAnn - SPARK_parent_proband_SNVs2$cnvSTARTAnn <= 10000, ] # only keep SNVs in CNVs <= 10,000 bp
SPARK_parent_proband_SNVs2$UID <- paste(SPARK_parent_proband_SNVs2$X.Sample, SPARK_parent_proband_SNVs2$X.id, sep='.')
SPARK_parent_proband_SNVs2 <- SPARK_parent_proband_SNVs2[!(SPARK_parent_proband_SNVs2$LoF & SPARK_parent_proband_SNVs2$freq_max > 0.01),] #remove LoF with snv freq > 1%

SPARK_parent_proband_ExonicSizes2 <- data.frame(SPARK_parent_proband_SNVsExonicSizes2[[2]])

# WGS3
SPARK_parent_proband_SNVsExonicSizes3 <- yaml::yaml.load_file("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/SPARKWGS3_parent_proband_SNVsExonicSizes.yaml")
SPARK_parent_proband_SNVs3 <- data.frame(SPARK_parent_proband_SNVsExonicSizes3[[1]])
SPARK_parent_proband_SNVs3 <- SPARK_parent_proband_SNVs3[SPARK_parent_proband_SNVs3$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')## filerting SNVs?
                                                         | SPARK_parent_proband_SNVs3$LoF,] 
SPARK_parent_proband_SNVs3 <- SPARK_parent_proband_SNVs23[SPARK_parent_proband_SNVs23$cnvENDAnn - SPARK_parent_proband_SNVs23$cnvSTARTAnn <= 10000, ] # only keep SNVs in CNVs <= 10,000 bp
SPARK_parent_proband_SNVs3$UID <- paste(SPARK_parent_proband_SNVs3$X.Sample, SPARK_parent_proband_SNVs3$X.id, sep='.')
SPARK_parent_proband_SNVs3 <- SPARK_parent_proband_SNVs3[!(SPARK_parent_proband_SNVs3$LoF & SPARK_parent_proband_SNVs3$freq_max > 0.01),] #remove LoF with snv freq > 1%

SPARK_parent_proband_ExonicSizes3 <- data.frame(SPARK_parent_proband_SNVsExonicSizes3[[2]])

# combine
SPARK_parent_proband_SNVs <- rbind(SPARK_parent_proband_SNVs1, SPARK_parent_proband_SNVs2, SPARK_parent_proband_SNVs3)
SPARK_parent_proband_ExonicSizes <- rbind(SPARK_parent_proband_ExonicSizes1, SPARK_parent_proband_ExonicSizes2, SPARK_parent_proband_ExonicSizes3)


## Get filtered metadata file:
metadata1 <- Filter_Metadata("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_1/SPARK_WGS_1_metadata_relfixed.tsv")
metadata2 <- Filter_Metadata("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_2/SPARK_WGS_2_metadata_relfixed.tsv")
metadata3 <- Filter_Metadata("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_3/SPARK_WGS_3_metadata_relfixed.tsv")
SPARK_meta <- rbind(metadata1, metadata2, metadata3) 
  

## Get SPARK IDs by relation
SPARK_father_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'father']
SPARK_mother_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'mother']
SPARK_male_proband_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'proband' &
                                                 SPARK_meta$Sex == "male"]
SPARK_female_proband_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'proband'&
                                                   SPARK_meta$Sex == "female"]

#####################################################################################################################################################################################################################################################

Get_CH_hit_By_Individual_ExonicSize <- function(SNVs, exonic_sizes,
                                                M.proband_IDs, F.proband_IDs,
                                                father_IDs, mother_IDs,
                                                cumulative = F,
                                                filter_SNV = T) {
  ## returns df with two cols: CH hit count and exonic size 
  if (filter_SNV) {
    SNVs <- SNVs[which(SNVs$gnomAD_pRec >= 0.9 & SNVs$gnomAD_oe_lof_upper >= 0.35),] 
  }
  
  df <- exonic_sizes %>% group_by(Sample.ID) # groups exonic sizes by sampleID
  if (!cumulative) {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% Sample.ID,]), # CH hits for each sample.ID
                           exSize=sum(exonicSize))
  }else {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% df$Sample.ID[df$cuts <= cuts[[1]]],])/
                             sum(df$exonicSize[df$cuts <= cuts[[1]]]))
  }
  
  df <- data.frame(df)
  df <- df[!is.nan(df$CH_hit),]
  
  ## add Relation column
  df$Relation <- NA
  df$Relation[df$Sample.ID %in% M.proband_IDs] <- 'Proband-male'
  df$Relation[df$Sample.ID %in% F.proband_IDs] <- 'Proband-female'
  
  df$Relation[df$Sample.ID %in% father_IDs] <- 'Father'
  df$Relation[df$Sample.ID %in% mother_IDs] <- 'Mother'
  df <- df[!is.na(df$Relation),] #######
  
  return (df)
}

#### SPARK 
SPARK_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(SNVs=SPARK_parent_proband_SNVs,
                                                            exonic_sizes=SPARK_parent_proband_ExonicSizes,
                                                            M.proband_IDs=SPARK_male_proband_IDs, 
                                                            F.proband_IDs=SPARK_female_proband_IDs,
                                                            SPARK_father_IDs, SPARK_mother_IDs,
                                                            filter_SNV = F)
SPARK_CH_hits_pRec <- Get_CH_hit_By_Individual_ExonicSize(SNVs=SPARK_parent_proband_SNVs,
                                                            exonic_sizes=SPARK_parent_proband_ExonicSizes,
                                                            M.proband_IDs=SPARK_male_proband_IDs, 
                                                            F.proband_IDs=SPARK_female_proband_IDs,
                                                            SPARK_father_IDs, SPARK_mother_IDs,
                                                            filter_SNV = T)


write.table(SPARK_CH_hits_nofilt, "./CH_count_del_size/data/CH_hits/SPARK_CH_hits_nofilt_sex.tsv",  
            sep="\t", row.names=F, quote=F, col.names=T)

write.table(SPARK_CH_hits_pRec, "./CH_count_del_size/data/CH_hits/SPARK_CH_hits_pRec_sex.tsv",  
            sep="\t", row.names=F, quote=F, col.names=T)

#####################################################################################################################################################################################################################
### Get plots #######

Get_CH_DelSize_Plots <- function(CH_hits_df, title){
  ## Saves no. CH events vs. exonic del size plot for CH_hits_df 
  
  controls.only <- CH_hits_df[which(CH_hits_df$Relation %in% c("Mother", "Father")),]
  CH_hits_df$Relation[which(!CH_hits_df$Relation %in% c("Mother", "Father"))] <- "Proband"
  
  ggplot(data=CH_hits_df, aes(x=exSize, y=CH_hit, group=1)) +
    geom_point(aes(color=Relation)) +
    geom_smooth(data=controls.only, method='lm', fullrange = T) +
    geom_abline(slope=sum(controls.only$CH_hit)/sum(controls.only$exSize)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    labs(y='No. of CH Events', x='Total deleted exonic size (bp)')+  
    ggtitle(sprintf("%s - No. CH Events v. Exonic Del Size", title)) 
  plot.path <- sprintf("./CH_count_del_size/figures/%s_CH_hits_smoothed.png", title)
  ggsave(plot.path, width = 7, height = 5)
}


##############################################################################################################################################################################################################################################
Get_CH_DelSize_Plots(SPARK_CH_hits_nofilt, "SPARK_nofilt")
Get_CH_DelSize_Plots(subset(SPARK_CH_hits_nofilt, exSize < 20000), "SPARK_nofilt_20kbless")

Get_CH_DelSize_Plots(SPARK_CH_hits_pRec, "SPARK_pRec")
Get_CH_DelSize_Plots(subset(SPARK_CH_hits_pRec, exSize < 20000), "SPARK_pRec_20kbless")
