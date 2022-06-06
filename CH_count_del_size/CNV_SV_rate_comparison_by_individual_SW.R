library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)

setwd("/Users/shaniawu/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/")

### Import & process data (from MSSNG+SSC_unbiased_burden_analysis.R) ### ##################################
Get_Failed_QC_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('./MSSNG+SSC/MSSNG+SSC_BA/data/excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) 
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('./MSSNG+SSC/MSSNG+SSC_BA/data/samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) 
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("./MSSNG+SSC/MSSNG+SSC_BA/data/MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("./MSSNG+SSC/MSSNG+SSC_BA/data/MSSNG+SSC.CNVs.tsv")
  ls_exclude <- unique(rbind(ls1, ls2))
  
  full_metadata <- fread(file, data.table = F)
  familyID_exclude <- full_metadata$`Family ID`[which(full_metadata$`Sample ID` %in% ls_exclude$Sample)]
  sampleID_exclude <- full_metadata$`Sample ID`[which(full_metadata$`Family ID` %in% familyID_exclude)]
  
  return (sampleID_exclude)
}

Get_Filtered_Metadata <- function(metadata_file, child_relation='proband', 
                                  comparison='parent-child') {
  meta <- read.delim(metadata_file, stringsAsFactors = F)
  
  # Get samples that failed QC
  failed_QC <- Get_Failed_QC_Samples()
  # Get samples with cr CNV
  samples_with_cr_CNV <- Get_cr_CNV_Samples()
  # Get families to exclude that have parents AS the proband as well
  families_to_exclude <- meta$Family.ID[(meta$Sample.ID[meta$Relation == 'proband'] %in% meta$Mother.ID) |
                                          (meta$Sample.ID[meta$Relation == 'proband'] %in% meta$Father.ID)]
  # Get samples with ASD associated CNVs and LOFs
  ASD_associated_samples <- Get_SampleID_ASD_Associated_CNVs_and_LOFs(metadata_file)
  
  # Filter out from metadata
  meta <- meta[!meta$Sample.ID %in% failed_QC, ]
  meta <- meta[!meta$Sample.ID %in% samples_with_cr_CNV, ]
  meta <- meta[!meta$Family.ID %in% families_to_exclude, ]
  meta <- meta[!meta$Sample.ID %in% ASD_associated_samples, ]
  
  child_relations <- ifelse(child_relation == 'proband', c('proband'), c('unaffected sibling', 'other sibling'))
  # Filter out irrelevant child if parent-child comparison
  if (comparison == 'parent-child') {
    meta <- meta[meta$Relation %in% child_relations | meta$Relation == 'mother'
                 | meta$Relation == 'father', ]
    # Get child-parent trios only
    IDs_relations <- meta[, c('Family.ID', 'Sample.ID', 'Relation')]
    IDs_relations$FamRel <- paste(IDs_relations$Family.ID, IDs_relations$Relation, sep='.')
    unique_famIDs_and_relations <- dplyr::distinct(IDs_relations,FamRel,
                                                   .keep_all=T)
    famID_freqs <- data.frame(table(unique_famIDs_and_relations$Family.ID))
    trios_IDs <- unique_famIDs_and_relations$Sample.ID[unique_famIDs_and_relations$Family.ID %in%
                                                         famID_freqs$Var1[famID_freqs$Freq == 3]] # 1 child, 1 father, 1 mother
    meta <- meta[meta$Sample.ID %in% trios_IDs, ]
  }
  # Filter out parents if child-child comparison
  else if (comparison == 'child-child') {
    meta <- meta[meta$Relation %in% c('proband', 'unaffected sibling', 'other sibling'), ]
    # Get only samples from families with 1 of each sibling
    proband_fam_IDs <- meta$Family.ID[meta$Relation == 'proband']
    sibling_fam_IDs <- meta$Family.ID[meta$Relation %in% c('unaffected sibling', 'other sibling')]
    meta <- meta[meta$Family.ID %in% proband_fam_IDs &
                   meta$Family.ID %in% sibling_fam_IDs,]
  }
  
  meta$Status <- ifelse(meta$Relation %in% child_relations, 1, 0)
  return (meta)
}

### MSSNG Parent-Proband 
MSSNG_metadata_path <- "./MSSNG+SSC/MSSNG+SSC_BA/data/MSSNG_metadata.tsv"
MSSNG_meta <- Get_Filtered_Metadata(MSSNG_metadata_path, child_relation = 'proband')

# MSSNG Proband, Mother, Father IDs
MSSNG_father_IDs <- MSSNG_meta$Sample.ID[MSSNG_meta$Relation == 'father']
MSSNG_mother_IDs <- MSSNG_meta$Sample.ID[MSSNG_meta$Relation == 'mother']
MSSNG_proband_IDs <- MSSNG_meta$Sample.ID[MSSNG_meta$Relation == 'proband']

# MSSNG CNV SNVS
MSSNG_parent_proband_SNVsExonicSizes <- yaml::yaml.load_file("./MSSNG+SSC/MSSNG+SSC_BA/data/MSSNG_parent_proband_SNVsExonicSizes.yaml")
MSSNG_parent_proband_SNVs <- data.frame(MSSNG_parent_proband_SNVsExonicSizes[[1]])
MSSNG_parent_proband_SNVs <- MSSNG_parent_proband_SNVs[MSSNG_parent_proband_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')## filerting SNVs?
                                                       | MSSNG_parent_proband_SNVs$LoF,] 
MSSNG_parent_proband_SNVs$UID <- paste(MSSNG_parent_proband_SNVs$X.Sample, MSSNG_parent_proband_SNVs$X.id, sep='.')
MSSNG_parent_proband_SNVs <- MSSNG_parent_proband_SNVs[!(MSSNG_parent_proband_SNVs$LoF & MSSNG_parent_proband_SNVs$freq_max > 0.01),] #remove LoF with snv freq > 1%
MSSNG_parent_proband_ExonicSizes <- data.frame(MSSNG_parent_proband_SNVsExonicSizes[[2]])

# MSSNG SV SNVS
MSSNG_parent_proband_SV_SNVsExonicSizes <- yaml::yaml.load_file("./MSSNG+SSC/MSSNG+SSC_BA/data/MSSNG_parent_proband_SV_SNVsExonicSizes.yaml")
MSSNG_parent_proband_SV_SNVs <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[1]])
MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[MSSNG_parent_proband_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')## filerting SNVs?
                                                             | MSSNG_parent_proband_SV_SNVs$LoF,] 
MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[MSSNG_parent_proband_SV_SNVs$cnvENDAnn - MSSNG_parent_proband_SV_SNVs$cnvSTARTAnn <= 10000, ] # only keep SNVs in CNVs <= 10,000 bp
MSSNG_parent_proband_SV_SNVs$UID <- paste(MSSNG_parent_proband_SV_SNVs$X.Sample, MSSNG_parent_proband_SV_SNVs$X.id, sep='.')
MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[!(MSSNG_parent_proband_SV_SNVs$LoF & MSSNG_parent_proband_SV_SNVs$freq_max > 0.01), ] #remove LoF with snv freq > 1%

MSSNG_parent_proband_SV_ExonicSizes <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[2]])

# Combine CNV and SV exonic sizes & SNVs
MSSNG_parent_proband_all_ExonicSizes <- merge(MSSNG_parent_proband_ExonicSizes, MSSNG_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
MSSNG_parent_proband_all_ExonicSizes$exonicSize.x[is.na(MSSNG_parent_proband_all_ExonicSizes$exonicSize.x)] <- 0
MSSNG_parent_proband_all_ExonicSizes$exonicSize.y[is.na(MSSNG_parent_proband_all_ExonicSizes$exonicSize.y)] <- 0
MSSNG_parent_proband_all_ExonicSizes$exonicSize <- pmax(MSSNG_parent_proband_all_ExonicSizes$exonicSize.x, MSSNG_parent_proband_all_ExonicSizes$exonicSize.y)
MSSNG_parent_proband_combined_SNVs <- rbind(MSSNG_parent_proband_SNVs, MSSNG_parent_proband_SV_SNVs)
MSSNG_parent_proband_proc_SNVs <- dplyr::distinct(MSSNG_parent_proband_combined_SNVs, UID, .keep_all=T)

# Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
MSSNG_parent_proband_SVs_ILMN <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[3]])
MSSNG_parent_proband_SVs_ILMN <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$length <= 10000
                                                               & MSSNG_parent_proband_SVs_ILMN$Sample.ID %in% MSSNG_meta$Sample.ID, ]
MSSNG_parent_proband_SVs_CG <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[4]])
MSSNG_parent_proband_SVs_CG <- MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$length <= 10000
                                                           & MSSNG_parent_proband_SVs_CG$Sample.ID %in% MSSNG_meta$Sample.ID, ]
MSSNG_homozyg_SVs <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
MSSNG_homozyg_SVs <- rbind(MSSNG_homozyg_SVs,MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')])
MSSNG_homozyg_SVs$SVUID <- with(MSSNG_homozyg_SVs, paste0(sample,CHROM,START,END))

MSSNG_parent_proband_proc_SNVs <- MSSNG_parent_proband_proc_SNVs[!with(MSSNG_parent_proband_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% MSSNG_homozyg_SVs$SVUID,]

MSSNG_parent_proband_proc_SNVs$Min_Dist <- by(MSSNG_parent_proband_proc_SNVs, seq_len(nrow(MSSNG_parent_proband_proc_SNVs)),
                                              function(r) r$MIN = min(abs(r$POS - MSSNG_parent_proband_proc_SNVs$POS[MSSNG_parent_proband_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                                       MSSNG_parent_proband_proc_SNVs$CHROM == r$CHROM &
                                                                                                                       MSSNG_parent_proband_proc_SNVs$POS != r$POS])) > 50)
MSSNG_parent_proband_proc_SNVs <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$Min_Dist,]

SNV_in_SV_only_MSSNG <- MSSNG_parent_proband_proc_SNVs[!MSSNG_parent_proband_proc_SNVs$UID %in% MSSNG_parent_proband_SNVs$UID, ]
SNV_in_CNV_only_MSSNG <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$UID %in% MSSNG_parent_proband_SNVs$UID, ]


### SSC Parent-Proband
SSC_metadata_path <- "./MSSNG+SSC/MSSNG+SSC_BA/data/SSC_metadata.tsv"
SSC_meta <- Get_Filtered_Metadata(SSC_metadata_path, child_relation = 'proband')

# SSC Proband, Father, Mother IDs
SSC_father_IDs <- SSC_meta$Sample.ID[SSC_meta$Relation == 'father']
SSC_mother_IDs <- SSC_meta$Sample.ID[SSC_meta$Relation == 'mother']
SSC_proband_IDs <- SSC_meta$Sample.ID[SSC_meta$Relation == 'proband']

SSC_parent_proband_SNVsExonicSizes <- yaml::yaml.load_file("./MSSNG+SSC/MSSNG+SSC_BA/data/SSC_parent_proband_SNVsExonicSizes.yaml")
SSC_parent_proband_SNVs <- data.frame(SSC_parent_proband_SNVsExonicSizes[[1]])
SSC_parent_proband_SNVs <- SSC_parent_proband_SNVs[SSC_parent_proband_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                   | SSC_parent_proband_SNVs$LoF,] 
SSC_parent_proband_SNVs$UID <- paste(SSC_parent_proband_SNVs$X.Sample, SSC_parent_proband_SNVs$X.id, sep='.')
SSC_parent_proband_SNVs <- SSC_parent_proband_SNVs[!(SSC_parent_proband_SNVs$LoF & SSC_parent_proband_SNVs$freq_max > 0.01), ]
SSC_parent_proband_ExonicSizes <- data.frame(SSC_parent_proband_SNVsExonicSizes[[2]])

SSC_parent_proband_SV_SNVsExonicSizes <- yaml::yaml.load_file("./MSSNG+SSC/MSSNG+SSC_BA/data/SSC_parent_proband_SV_SNVsExonicSizes.yaml")
SSC_parent_proband_SV_SNVs <- data.frame(SSC_parent_proband_SV_SNVsExonicSizes[[1]])
SSC_parent_proband_SV_SNVs <- SSC_parent_proband_SV_SNVs[SSC_parent_proband_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                         | SSC_parent_proband_SV_SNVs$LoF,] 
SSC_parent_proband_SV_SNVs <- SSC_parent_proband_SV_SNVs[SSC_parent_proband_SV_SNVs$cnvENDAnn - SSC_parent_proband_SV_SNVs$cnvSTARTAnn <= 10000, ]
SSC_parent_proband_SV_SNVs$UID <- paste(SSC_parent_proband_SV_SNVs$X.Sample, SSC_parent_proband_SV_SNVs$X.id, sep='.')
SSC_parent_proband_SV_SNVs <- SSC_parent_proband_SV_SNVs[!(SSC_parent_proband_SV_SNVs$LoF & SSC_parent_proband_SV_SNVs$freq_max > 0.01), ]
SSC_parent_proband_SV_ExonicSizes <- data.frame(SSC_parent_proband_SV_SNVsExonicSizes[[2]])

SSC_parent_proband_all_ExonicSizes <- merge(SSC_parent_proband_ExonicSizes, SSC_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T)
# Combine 
SSC_parent_proband_all_ExonicSizes <- merge(SSC_parent_proband_ExonicSizes, SSC_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
SSC_parent_proband_all_ExonicSizes$exonicSize.x[is.na(SSC_parent_proband_all_ExonicSizes$exonicSize.x)] <- 0
SSC_parent_proband_all_ExonicSizes$exonicSize.y[is.na(SSC_parent_proband_all_ExonicSizes$exonicSize.y)] <- 0
SSC_parent_proband_all_ExonicSizes$exonicSize <- pmax(SSC_parent_proband_all_ExonicSizes$exonicSize.x, SSC_parent_proband_all_ExonicSizes$exonicSize.y)
SSC_parent_proband_combined_SNVs <- rbind(SSC_parent_proband_SNVs, SSC_parent_proband_SV_SNVs)
SSC_parent_proband_proc_SNVs <- dplyr::distinct(SSC_parent_proband_combined_SNVs, UID, .keep_all=T)

ssc_all_SNVs_combined <- rbind(SSC_parent_proband_proc_SNVs, SSC_parent_US_proc_SNVs)
ssc_all_SNVs_combined_unique <- dplyr::distinct(ssc_all_SNVs_combined, UID, .keep_all = T)
ssc_all_ExonicSizes_combined <- rbind(SSC_parent_proband_all_ExonicSizes, SSC_parent_US_all_ExonicSizes)
ssc_all_ExonicSizes_combined_unique <- dplyr::distinct(ssc_all_ExonicSizes_combined, Sample.ID, .keep_all = T)
SNV_in_SV_only_SSC <- ssc_all_SNVs_combined_unique[!ssc_all_SNVs_combined_unique$UID %in% SSC_parent_proband_SNVs$UID &
                                                     !ssc_all_SNVs_combined_unique$UID %in% SSC_parent_US_SNVs$UID, ]

# Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
SSC_parent_proband_SVs <- data.frame(SSC_parent_proband_SV_SNVsExonicSizes[[3]])
SSC_parent_proband_SVs <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$length <= 10000
                                                 & SSC_parent_proband_SVs$Sample.ID %in% SSC_meta$Sample.ID, ]

SSC_homozyg_SVs <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
SSC_homozyg_SVs$SVUID <- with(SSC_homozyg_SVs, paste0(sample,CHROM,START,END))

SSC_parent_proband_proc_SNVs <- SSC_parent_proband_proc_SNVs[!with(SSC_parent_proband_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% SSC_homozyg_SVs$SVUID,]

SSC_parent_proband_proc_SNVs$Min_Dist <- by(SSC_parent_proband_proc_SNVs, seq_len(nrow(SSC_parent_proband_proc_SNVs)),
                                            function(r) r$MIN = min(abs(r$POS - SSC_parent_proband_proc_SNVs$POS[SSC_parent_proband_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                                   SSC_parent_proband_proc_SNVs$CHROM == r$CHROM &
                                                                                                                   SSC_parent_proband_proc_SNVs$POS != r$POS])) > 50)
SSC_parent_proband_proc_SNVs <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$Min_Dist,]


#####################################################################################################################################################################################################################################################

Get_CH_hit_By_Individual_ExonicSize <- function(SNVs,
                                                exonic_sizes,
                                                proband_IDs,
                                                father_IDs,
                                                mother_IDs,
                                                cumulative = F,
                                                filter_SNV = T) {
  ## returns df with two cols: CH hit count and exonic size 
  if (filter_SNV) {
    SNVs <- SNVs[which(SNVs$gnomAD_pRec >= 0.9 & SNVs$gnomAD_oe_lof_upper >= 0.35),] # changed from Faraz
  }

  df <- exonic_sizes %>% group_by(Sample.ID) # groups exonic sizes by sampleID
  if (!cumulative) {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% Sample.ID,]), # CH hits for each sample.ID
                           exSize=sum(exonicSize))
  }
  else {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% df$Sample.ID[df$cuts <= cuts[[1]]],])/
                             sum(df$exonicSize[df$cuts <= cuts[[1]]]))
  }
  
  df <- data.frame(df)
  df <- df[!is.nan(df$CH_hit),]
  
  ## add Relation column
  df$Relation <- NA
  df$Relation[df$Sample.ID %in% proband_IDs] <- 'Proband'
  df$Relation[df$Sample.ID %in% father_IDs] <- 'Father'
  df$Relation[df$Sample.ID %in% mother_IDs] <- 'Mother'
  df <- df[!is.na(df$Relation),]
  
  return (df)
}

#####################################################################################################################################################################################################################
## Get plots 

Get_CH_DelSize_Plots <- function(CH_hits_df, title){
  ## Saves no. CH events vs. exonic del size plot for CH_hits_df 
  
  ggplot(data=CH_hits_df, aes(x=exSize, y=CH_hit, group=1)) +
    geom_point(aes(color=Relation)) +
    geom_smooth(method='lm') +
    geom_abline(slope=sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    labs(y='No. of CH Events', x='Total deleted exonic size (bp)') +
    ggtitle(sprintf("%s- No. CH Events v. Exonic Del Size", title)) 
  plot.path <- sprintf("./CH_count_del_size/figures/%s_CH_hits_nofilt_smoothed.png", title)
  ggsave(plot.path, width = 7, height = 5)
}


#### MSSNG ####
MSSNG_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs,
                                                            MSSNG_parent_proband_all_ExonicSizes,
                                                            MSSNG_proband_IDs, MSSNG_father_IDs, MSSNG_mother_IDs,
                                                            filter_SNV = F)
Get_CH_DelSize_Plots(MSSNG_CH_hits_nofilt, "MSSNG")

MSSNG_CH_hits_nofilt_less75kb <- MSSNG_CH_hits_nofilt[(which(MSSNG_CH_hits_nofilt$exSize <= 75000)),] # restrict del size to < 75 kb
Get_CH_DelSize_Plots(MSSNG_CH_hits_nofilt_less75kb, "MSSNG_<75kb")

## Plot by variant type (effect_priority) - MSSNG
# MSSNG nonsynonymous
nonsyn.MSSNG.SNVs <- MSSNG_parent_proband_proc_SNVs[which(MSSNG_parent_proband_proc_SNVs$effect_priority == "nonsynonymous SNV"),]
nonsyn.MSSNG_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(nonsyn.MSSNG.SNVs, MSSNG_parent_proband_all_ExonicSizes, 
                                                            MSSNG_proband_IDs, MSSNG_father_IDs, MSSNG_mother_IDs, 
                                                            filter_SNV = F)
Get_CH_DelSize_Plots(nonsyn.MSSNG_CH_hits, "MSSNG_nonsyn")

nonsyn.MSSNG_CH_hit_less75kb <- nonsyn.MSSNG_CH_hits[(which(nonsyn.MSSNG_CH_hits$exSize <= 75000)),] # restrict del size to < 75 kb
Get_CH_DelSize_Plots(nonsyn.MSSNG_CH_hit_less75kb, "MSSNG_nonsyn_<75kb")

# MSSNG synonymous
syn.MSSNG.SNVs <- MSSNG_parent_proband_proc_SNVs[which(MSSNG_parent_proband_proc_SNVs$effect_priority == "synonymous SNV"),]
syn.MSSNG_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(syn.MSSNG.SNVs, MSSNG_parent_proband_all_ExonicSizes, 
                                                            MSSNG_proband_IDs, MSSNG_father_IDs, MSSNG_mother_IDs, 
                                                            filter_SNV = F)
Get_CH_DelSize_Plots(syn.MSSNG_CH_hits, "MSSNG_syn")

syn.MSSNG_CH_hit_less75kb <- syn.MSSNG_CH_hits[(which(syn.MSSNG_CH_hits$exSize <= 75000)),] # restrict del size to < 75 kb
Get_CH_DelSize_Plots(syn.MSSNG_CH_hit_less75kb, "MSSNG_syn_<75kb")

#### SSC ####
SSC_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs, 
                                                          SSC_parent_proband_all_ExonicSizes,
                                                          SSC_proband_IDs, SSC_father_IDs, SSC_mother_IDs,
                                                          filter_SNV = F)
Get_CH_DelSize_Plots(SSC_CH_hits_nofilt, "SSC")

SSC_CH_hits_nofilt_less25kb <- SSC_CH_hits_nofilt[(which(SSC_CH_hits_nofilt$exSize <= 25000)),] # restrict del size to < 25 kb
Get_CH_DelSize_Plots(SSC_CH_hits_nofilt_less25kb, "SSC_<25kb")


## Plot by variant type (effect_priority) - SSC

# SSC nonsynonymous
nonsyn.SSC.SNVs <- SSC_parent_proband_proc_SNVs[which(SSC_parent_proband_proc_SNVs$effect_priority == "nonsynonymous SNV"),]
nonsyn.SSC_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(nonsyn.SSC.SNVs, SSC_parent_proband_all_ExonicSizes, 
                                                            SSC_proband_IDs, SSC_father_IDs, SSC_mother_IDs, 
                                                            filter_SNV = F)
Get_CH_DelSize_Plots(nonsyn.SSC_CH_hits, "SSC_nonsyn")

nonsyn.SSC_CH_hits_less25kb <- nonsyn.SSC_CH_hits[(which(nonsyn.SSC_CH_hits$exSize <= 25000)),] # restrict del size to < 25 kb
Get_CH_DelSize_Plots(nonsyn.SSC_CH_hits_less25kb, "SSC_nonsyn_<25kb")

# SSC synonymous
syn.SSC.SNVs <- SSC_parent_proband_proc_SNVs[which(SSC_parent_proband_proc_SNVs$effect_priority == "synonymous SNV"),]
syn.SSC_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(syn.SSC.SNVs, SSC_parent_proband_all_ExonicSizes, 
                                                         SSC_proband_IDs, SSC_father_IDs, SSC_mother_IDs, 
                                                         filter_SNV = F)
Get_CH_DelSize_Plots(syn.SSC_CH_hits, "SSC_syn")

syn.SSC_CH_hits_less25kb <- syn.SSC_CH_hits[(which(syn.SSC_CH_hits$exSize <= 25000)),] # restrict del size to < 25 kb
Get_CH_DelSize_Plots(syn.SSC_CH_hits_less25kb, "SSC_syn_<25kb")


#### Combined MSSNG and SSC datasets ####

## No pRec filter:####
MSSNG_SSC_CH_hits_nofilt <- rbind(MSSNG_CH_hits_nofilt, SSC_CH_hits_nofilt)

Get_CH_DelSize_Plots(MSSNG_SSC_CH_hits_nofilt, "MSSNG+SSC")


MSSNG_SSC_CH_hits_nofilt_less75kb <- MSSNG_SSC_CH_hits_nofilt[(which(MSSNG_SSC_CH_hits_nofilt$exSize <= 75000)),] # restrict del size to < 75 kb
Get_CH_DelSize_Plots(MSSNG_SSC_CH_hits_nofilt_less75kb, "MSSNG+SSC_<75kb")


## Plot by variant type (effect_priority) - MSSNG+SSC
MSSNG_SSC_parent_proband_proc_SNVs <- rbind(MSSNG_parent_proband_proc_SNVs, SSC_parent_proband_proc_SNVs)
MSSNG_SSC_parent_proband_all_ExonicSize <- rbind(MSSNG_parent_proband_all_ExonicSizes, SSC_parent_proband_all_ExonicSizes)

# MSSNG+SSC nonsynonymous
nonsyn.SNVs <- MSSNG_SSC_parent_proband_proc_SNVs[which(MSSNG_SSC_parent_proband_proc_SNVs$effect_priority == "nonsynonymous SNV"),]
nonsyn.MSSNG_SSC_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(nonsyn.SNVs, MSSNG_SSC_parent_proband_all_ExonicSize, 
                                                                c(MSSNG_proband_IDs, SSC_proband_IDs), 
                                                                c(MSSNG_father_IDs, SSC_father_IDs), 
                                                                c(MSSNG_mother_IDs, SSC_mother_IDs),
                                                                filter_SNV = F)
Get_CH_DelSize_Plots(nonsyn.MSSNG_SSC_CH_hits, "MSSNG+SSC_nonsyn")

nonsyn.MSSNG_SSC_CH_hits_less75kb <- nonsyn.MSSNG_SSC_CH_hits[(which(nonsyn.MSSNG_SSC_CH_hits$exSize <= 75000)),] # restrict del size to < 75 kb
Get_CH_DelSize_Plots(nonsyn.MSSNG_SSC_CH_hits_less75kb, "MSSNG+SSC_nonsyn_<75kb")

# MSSNG+SSC synonymous
syn.SNVs <- MSSNG_SSC_parent_proband_proc_SNVs[which(MSSNG_SSC_parent_proband_proc_SNVs$effect_priority == "synonymous SNV"),]
syn.MSSNG_SSC_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(syn.SNVs, MSSNG_SSC_parent_proband_all_ExonicSize, 
                                                                c(MSSNG_proband_IDs, SSC_proband_IDs), 
                                                                c(MSSNG_father_IDs, SSC_father_IDs), 
                                                                c(MSSNG_mother_IDs, SSC_mother_IDs),
                                                                filter_SNV = F)
Get_CH_DelSize_Plots(syn.MSSNG_SSC_CH_hits, "MSSNG+SSC_syn")

syn.MSSNG_SSC_CH_hits_less75kb <- syn.MSSNG_SSC_CH_hits[(which(syn.MSSNG_SSC_CH_hits$exSize <= 75000)),] # restrict del size to < 75 kb
Get_CH_DelSize_Plots(syn.MSSNG_SSC_CH_hits_less75kb, "MSSNG+SSC_syn_<75kb")


## MSSNG+SSC combined with pRec > 0.9 filter ####
MSSNG_CH_hits_pRec <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs,
                                                          MSSNG_parent_proband_all_ExonicSizes,
                                                          MSSNG_proband_IDs, MSSNG_father_IDs, MSSNG_mother_IDs,
                                                          filter_SNV = T) 
SSC_CH_hits_pRec <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs,
                                                        SSC_parent_proband_all_ExonicSizes,
                                                        SSC_proband_IDs, SSC_father_IDs, SSC_mother_IDs,
                                                        filter_SNV = T)
MSSNG_SSC_CH_hits_pRec <- rbind(MSSNG_CH_hits_pRec, SSC_CH_hits_pRec)
Get_CH_DelSize_Plots(MSSNG_SSC_CH_hits_pRec, "MSSNG+SSC_pRec0.9")

# MSSNG_SSC_CH_hits_pRec_less75kb <- MSSNG_SSC_CH_hits_pRec[(which(MSSNG_SSC_CH_hits_pRec$exSize <= 75000)),] # restrict del size to < 75 kb
# Get_CH_DelSize_Plots(MSSNG_SSC_CH_hits_pRec_less75kb, "MSSNG+SSC_pRec0.9_<75kb")



##############################################################################################################################################################################################################################################
#### Fisher's Exact Tests ####

Get_Fisher_Res <- function(CH_hits_df, CH_hits_df_type, size_bin = F, pRec_size_bins = F) {
  ## Returns Fisher's exact test results as a dataframe, given a dataframe of CH hits and exonic sizes
  ## CH_hits_df: all CH hits 
  ## CH_hits_df_type: target type (nonsyn, syn) of CH hits
  
  ## All del sizes
  if (size_bin == F & pRec_size_bins == F){ 
    ## probands
    CH_hits_probands_target <- CH_hits_df_type[which(CH_hits_df_type$Relation == 'Proband'),]
    ## controls (father & mother)
    CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother'),]
    
    lm.model <- lm(CH_hit ~ exSize, CH_hits_df_type)
    intercept <- lm.model$coefficients[1]
    slope <- lm.model$coefficients[2]
    proband_resi <- CH_hits_probands_target$CH_hit - ((slope * CH_hits_probands_target$exSize) + intercept)
    control_resi <- CH_hits_controls_target$CH_hit - ((slope * CH_hits_controls_target$exSize) + intercept)
    # proband_resi <- CH_hits_probands_target$CH_hit - ((sum(CH_hits_df_type$CH_hit)/sum(CH_hits_df$exSize)) * CH_hits_probands_target$exSize)
    # control_resi <- CH_hits_controls_target$CH_hit - ((sum(CH_hits_df_type$CH_hit)/sum(CH_hits_df$exSize)) * CH_hits_controls_target$exSize)
  }
  
  ## Different exonic deletion size bins (non-pRec size bins)
  if (!size_bin == F){
    if (size_bin == "0-1kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & CH_hits_df_type$exSize <= 1000,]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & CH_hits_df_type$exSize <= 1000,]
    }
    if (size_bin == "1-2kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                              (CH_hits_df_type$exSize > 1000 & CH_hits_df_type$exSize <= 2000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                              (CH_hits_df_type$exSize > 1000 & CH_hits_df_type$exSize <= 2000),]
    }
    if (size_bin == "2-5kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                              (CH_hits_df_type$exSize > 2000 & CH_hits_df_type$exSize <= 5000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                              (CH_hits_df_type$exSize > 2000 & CH_hits_df_type$exSize <= 5000),]
    }
    if (size_bin == "5-10kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                              (CH_hits_df_type$exSize > 5000 & CH_hits_df_type$exSize <= 10000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                              (CH_hits_df_type$exSize > 5000 & CH_hits_df_type$exSize <= 10000),]
    }
    if (size_bin == "10kb+"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & CH_hits_df_type$exSize > 10000,]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & CH_hits_df_type$exSize > 10000,]
    }
    lm.model <- lm(CH_hit ~ exSize, CH_hits_df_type)
    intercept <- lm.model$coefficients[1]
    slope <- lm.model$coefficients[2]
    proband_resi <- CH_hits_probands_target$CH_hit - ((slope * CH_hits_probands_target$exSize) + intercept)
    control_resi <- CH_hits_controls_target$CH_hit - ((slope * CH_hits_controls_target$exSize) + intercept)
    # proband_resi <- CH_hits_probands_target$CH_hit - ((sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) * CH_hits_probands_target$exSize)
    # control_resi <- CH_hits_controls_target$CH_hit - ((sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) * CH_hits_controls_target$exSize)
  }
  ## Different exonic deletion size bins (pRec size bins)
  if (!pRec_size_bins == F){
    if (pRec_size_bins == "0-1kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & CH_hits_df_type$exSize <= 1000,]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & CH_hits_df_type$exSize <= 1000,]
    }
    if (pRec_size_bins == "1-5kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                                   (CH_hits_df_type$exSize > 1000 & CH_hits_df_type$exSize <= 5000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                                   (CH_hits_df_type$exSize > 1000 & CH_hits_df_type$exSize <= 5000),]
    }
    if (pRec_size_bins == "5-10kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                                   (CH_hits_df_type$exSize > 5000 & CH_hits_df_type$exSize <= 10000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                                   (CH_hits_df_type$exSize > 5000 & CH_hits_df_type$exSize <= 10000),]
    }
    if (pRec_size_bins == "10-20kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                                   (CH_hits_df_type$exSize > 10000 & CH_hits_df_type$exSize <= 20000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                                   (CH_hits_df_type$exSize > 10000 & CH_hits_df_type$exSize <= 20000),]
    }
    if (pRec_size_bins == "20-50kb"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & 
                                                   (CH_hits_df_type$exSize > 20000 & CH_hits_df_type$exSize <= 50000),]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & 
                                                   (CH_hits_df_type$exSize > 20000 & CH_hits_df_type$exSize <= 50000),]
    }
    if (pRec_size_bins == "50kb+"){
      CH_hits_probands_target <- CH_hits_df_type[CH_hits_df_type$Relation == 'Proband' & CH_hits_df_type$exSize > 50000,]
      CH_hits_controls_target <- CH_hits_df_type[CH_hits_df_type$Relation %in% c('Father','Mother') & CH_hits_df_type$exSize > 50000,]
    }
    lm.model <- lm(CH_hit ~ exSize, CH_hits_df_type)
    intercept <- lm.model$coefficients[1]
    slope <- lm.model$coefficients[2]
    proband_resi <- CH_hits_probands_target$CH_hit - ((slope * CH_hits_probands_target$exSize) + intercept)
    control_resi <- CH_hits_controls_target$CH_hit - ((slope * CH_hits_controls_target$exSize) + intercept)
    # proband_resi <- CH_hits_probands_target$CH_hit - ((sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) * CH_hits_probands_target$exSize)
    # control_resi <- CH_hits_controls_target$CH_hit - ((sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) * CH_hits_controls_target$exSize)
  }
    
  ## Make Fisher Exact Test:
  proband_resi <- proband_resi[!is.na(proband_resi)]
  proband_above_expect <- sum(proband_resi > 0)
  proband_below_expect <- sum(proband_resi < 0)
  
  control_resi <- control_resi[!is.na(control_resi)]
  control_above_expect <- sum(control_resi > 0)
  control_below_expect <- sum(control_resi < 0)
  
  ratio_case_control_above_expect <- proband_above_expect/control_above_expect
  ratio_case_control_below_expect <- proband_below_expect/control_below_expect
  case_control_expect_OR <- ratio_case_control_above_expect/ratio_case_control_below_expect
  case_control_fisher_df <- data.frame(Above=c(proband_above_expect,control_above_expect),
                                       Below=c(proband_below_expect,control_below_expect))
  rownames(case_control_fisher_df) <- c('Case','Control')
  case_control_fisher_res <- fisher.test(case_control_fisher_df, alternative='greater')
  case_control_fisher_res.df <- broom::tidy(case_control_fisher_res)
  
  ## add counts to fisher results
  case_control_fisher_res.df$case.above <- proband_above_expect
  case_control_fisher_res.df$case.below <- proband_below_expect
  case_control_fisher_res.df$control.above <- control_above_expect
  case_control_fisher_res.df$control.below <- control_below_expect
  
  return(case_control_fisher_res.df)
  
  file.path <- sprintf("./CH_count_del_size/data/%s.%s.fisher.tsv", CH_hits_df_type, pRec_size_bins)
  # write.csv(case_control_fisher_res.df, file.path)
}

#### Fisher's Test for MSSNG (all variants) #### 
MSSNG_all_fisher <- Get_Fisher_Res(MSSNG_CH_hits_nofilt, MSSNG_CH_hits_nofilt)
write.table(MSSNG_all_fisher, "./CH_count_del_size/data/MSSNG_all_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)

#### Fisher's Test for SSC (all variants) ####
SSC_all_fisher <- Get_Fisher_Res(SSC_CH_hits_nofilt, SSC_CH_hits_nofilt)
write.table(SSC_all_fisher, "./CH_count_del_size/data/SSC_all_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)

#### Fisher's Test for MSSNG+SSC (all variants) ####
MSSNG_SSC_all_fisher <- Get_Fisher_Res(MSSNG_SSC_CH_hits_nofilt, MSSNG_SSC_CH_hits_nofilt)
write.table(MSSNG_SSC_all_fisher, "./CH_count_del_size/data/MSSNG_SSC_all_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)

## MSSNG+SSC nonsyn
MSSNG_SSC_nonsyn_fisher <- Get_Fisher_Res(MSSNG_SSC_CH_hits_nofilt, nonsyn.MSSNG_SSC_CH_hits)
write.table(MSSNG_SSC_nonsyn_fisher, "./CH_count_del_size/data/MSSNG_SSC_nonsyn_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)


## MSSNG+SSC syn
MSSNG_SSC_syn_fisher <- Get_Fisher_Res(MSSNG_SSC_CH_hits_nofilt, syn.MSSNG_SSC_CH_hits)
write.table(MSSNG_SSC_syn_fisher, "./CH_count_del_size/data/MSSNG_SSC_syn_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)

#### MSSNG+SSC pRec > 0.9
# MSSNG_SSC_pRec_fisher <- Get_Fisher_Res(MSSNG_SSC_CH_hits_pRec, MSSNG_SSC_CH_hits_pRec, pRec_size_bins = T)
# write.table(MSSNG_SSC_all_fisher, "./CH_count_del_size/data/MSSNG_SSC_all_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)


#### Fisher's Test for Exonic Del Size Bins ####
  
Get_Fisher_Res_SizeBins <- function(CH_hits_df, name, pRec_size_bins = F){
  ## Writes a table of Fishers test results for 5 exonic deletion size bins under name
  fisher.res.comb <- data.frame()
  if (pRec_size_bins == F){
    for (size_bin in c("0-1kb", "1-2kb", "2-5kb", "5-10kb", "10kb+")){
      fisher.res <- Get_Fisher_Res(CH_hits_df, CH_hits_df, size_bin = size_bin)
      fisher.res <- cbind(size.bin = size_bin, fisher.res)
      fisher.res.comb <- rbind(fisher.res.comb, fisher.res)
      
      file.path <- sprintf("./CH_count_del_size/data/%s_sizes_fisher.tsv", name)
      write.table(fisher.res.comb, file.path, sep="\t", row.names=F, quote=F, col.names=T)
    }
  }
  if (pRec_size_bins == T){
    for (size_bin in c("0-1kb", "1-5kb", "5-10kb", "10-20kb", "20-50kb", "50kb+")){
      fisher.res <- Get_Fisher_Res(CH_hits_df, CH_hits_df, pRec_size_bins = size_bin)
      fisher.res <- cbind(size.bin = size_bin, fisher.res)
      fisher.res.comb <- rbind(fisher.res.comb, fisher.res)
      
      file.path <- sprintf("./CH_count_del_size/data/%s_sizes_fisher.tsv", name)
      write.table(fisher.res.comb, file.path, sep="\t", row.names=F, quote=F, col.names=T)
    }
  }
}

## MSSNG + SSC
MSSNG_SSC_all_fisher_sizes <- Get_Fisher_Res_SizeBins(MSSNG_SSC_CH_hits_nofilt, "MSSNG.SSC.all.nofilt")
MSSNG_SSC_nonsyn_fisher_sizes <- Get_Fisher_Res_SizeBins(nonsyn.MSSNG_SSC_CH_hits, "MSSNG.SSC.nonsyn.nofilt")
MSSNG_SSC_syn_fisher_sizes <- Get_Fisher_Res_SizeBins(syn.MSSNG_SSC_CH_hits, "MSSNG.SSC.syn.nofilt")
# pRec > 0.9
MSSNG_SSC_pRec_fisher_sizes <- Get_Fisher_Res_SizeBins(MSSNG_SSC_CH_hits_pRec, "MSSNG.SSC.pRec", pRec_size_bins = T) 

## MSSNG 
MSSNG_all_fisher_sizes <- Get_Fisher_Res_SizeBins(MSSNG_CH_hits_nofilt, "MSSNG.all.nofilt")
MSSNG_nonsyn_fisher_sizes <- Get_Fisher_Res_SizeBins(nonsyn.MSSNG_CH_hits, "MSSNG.nonsyn.nofilt")
MSSNG_syn_fisher_sizes <- Get_Fisher_Res_SizeBins(syn.MSSNG_CH_hits, "MSSNG.syn.nofilt")

## SSC 
SSC_all_fisher_sizes <- Get_Fisher_Res_SizeBins(SSC_CH_hits_nofilt, "SSC.all.nofilt")
MSSNG_SSC_nonsyn_fisher_sizes <- Get_Fisher_Res_SizeBins(nonsyn.SSC_CH_hits, "SSC.nonsyn.nofilt")
MSSNG_SSC_syn_fisher_sizes <- Get_Fisher_Res_SizeBins(syn.SSC_CH_hits, "SSC.syn.nofilt")








