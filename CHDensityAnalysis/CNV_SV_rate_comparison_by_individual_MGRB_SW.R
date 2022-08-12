################################################################################################################################################################################
# 
# CNV_SV_rate_comparison_by_individual_MGRB_SW.R
# purpose: - creates CH_hits tables (CH hit and total deleted exonic size tables)
#          - plots MGRB linear del size vs. CH count plots
#          - outputs MGRB Fisher's test results using linear observed line as cut-off 
# input: MGRB metadta, CRVs, CNV_SNVsExonicSizes_MGRB.yaml
# output: MGRB del size vs. CH count plots & Fisher's test results using linear 
#           del size vs. CH count as cut-off. e.g. CH_hits_nofilt_noCRVs_smoothed.png
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

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/")

### Import & process data  ######################################################################################################

# MGRB CRVs (to exclude)
CRVs <- readLines("./MGRB/data_processing/data/samples.with.crCNV.MGRB.txt")

# MGRB metadata
MGRB_metadata_path <- "~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MGRB/MGRB_set1_metadata.tsv"
MGRB_meta <- read.delim(MGRB_metadata_path, stringsAsFactors = F)
MGRB_meta <- MGRB_meta[!MGRB_meta$Sample.ID %in% CRVs,]

MGRB.F.sampleIDs <- MGRB_meta$Sample.ID[which(MGRB_meta$Sex == "female")]
MGRB.M.sampleIDs <- MGRB_meta$Sample.ID[which(MGRB_meta$Sex == "male")]

# MGRB SNVS
MGRB_SNVsExonicSizes <- yaml::yaml.load_file("./MGRB/data_processing/data/CNV_SNVsExonicSizes_MGRB.yaml")
MGRB_SNVs <- data.frame(MGRB_SNVsExonicSizes[[1]])
MGRB_SNVs <- MGRB_SNVs[MGRB_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV') | MGRB_SNVs$LoF,] 
MGRB_SNVs$UID <- paste(MGRB_SNVs$X.Sample, MGRB_SNVs$Original_VCFKEY, sep='.') 
MGRB_SNVs <- MGRB_SNVs[!(MGRB_SNVs$LoF & MGRB_SNVs$freq_max > 0.01),] # remove LoF with snv freq > 1%

MGRB_ExonicSizes <- data.frame(MGRB_SNVsExonicSizes[[2]])
MGRB_ExonicSizes.F <- MGRB_ExonicSizes[which(MGRB_ExonicSizes$Sample.ID %in% MGRB.F.sampleIDs),]
MGRB_ExonicSizes.M <- MGRB_ExonicSizes[which(MGRB_ExonicSizes$Sample.ID %in% MGRB.M.sampleIDs),]

MGRB_proc_SNVs <- dplyr::distinct(MGRB_SNVs, UID, .keep_all=T) 

# Filter out homozygous dels or FP SNV calls (< 50bp distance SNV)
MGRB_CNVs <- data.frame(MGRB_SNVsExonicSizes[[3]])
MGRB_CNVs <- MGRB_CNVs[MGRB_CNVs$length <= 10000 & MGRB_CNVs$Sample.ID %in% MGRB_meta$Sample.ID, ]
MGRB_CNVs <- MGRB_CNVs[which(MGRB_CNVs$SAMPLE.CN == 1),] # filter out homozygous deletions

MGRB_proc_SNVs$Min_Dist <- by(MGRB_proc_SNVs, seq_len(nrow(MGRB_proc_SNVs)),
                                              function(r) r$MIN = min(abs(r$POS - MGRB_proc_SNVs$POS[MGRB_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                                       MGRB_proc_SNVs$CHROM == r$CHROM &
                                                                                                                       MGRB_proc_SNVs$POS != r$POS])) > 50)
MGRB_proc_SNVs <- MGRB_proc_SNVs[MGRB_proc_SNVs$Min_Dist,]

# split into sexes
MGRB_proc_SNVs_F <- MGRB_proc_SNVs[which(MGRB_proc_SNVs$X.Sample %in% MGRB.F.sampleIDs),]
MGRB_proc_SNVs_M <- MGRB_proc_SNVs[which(MGRB_proc_SNVs$X.Sample %in% MGRB.M.sampleIDs),]


#####################################################################################################################################################################################################################################################

Get_CH_hit_By_Individual_ExonicSize <- function(SNVs,
                                                exonic_sizes,
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
  }
  else {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% df$Sample.ID[df$cuts <= cuts[[1]]],])/
                             sum(df$exonicSize[df$cuts <= cuts[[1]]]))
  }
  df <- data.frame(df)
  df <- df[!is.nan(df$CH_hit),]
  
  return (df)
}

#####################################################################################################################################################################################################################
## Get plots 

Get_CH_DelSize_Plots <- function(CH_hits_df, title){
  ## Saves no. CH events vs. exonic del size plot for CH_hits_df 
  ggplot(data=CH_hits_df, aes(x=exSize, y=CH_hit, group=1)) +
    geom_point() +
    geom_smooth(method='lm') +
    geom_abline(slope=sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    labs(y='No. of CH Events', x='Total deleted exonic size (bp)') +
    ggtitle(sprintf("%s- No. CH Events v. Exonic Del Size", title))
  plot.path <- sprintf("./CH_count_del_size/figures/%s_CH_hits_nofilt_noCRVs_smoothed.png", title)
  ggsave(plot.path, width = 7, height = 5)
}

Get_CH_DelSize_Sex_Plots <- function(CH_hits_df, CH_hits_df.F, CH_hits_df.M, title){
  ## Saves no. CH events vs. exonic del size plot for CH_hits_df 
  ggplot(data=CH_hits_df, aes(x=exSize, y=CH_hit, group=1)) +
    geom_point() +
    geom_smooth(data = CH_hits_df.F, method='lm', colour = "#FF9999") + # female observed line
    geom_smooth(data = CH_hits_df.M, method='lm', colour = "#56B4E9") + # male observed line
    geom_abline(slope=sum(CH_hits_df$CH_hit)/sum(CH_hits_df$exSize)) + # expected line
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    labs(y='No. of CH Events', x='Total deleted exonic size (bp)') +
    ggtitle(sprintf("%s- No. CH Events v. Exonic Del Size", title)) 
  plot.path <- sprintf("./CH_count_del_size/figures/%s_CH_hits_nofilt_noCRVs_smoothed.png", title)
  ggsave(plot.path, width = 7, height = 5)
}

#### MGRB ####
MGRB_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(MGRB_proc_SNVs,
                                                           MGRB_ExonicSizes,
                                                           filter_SNV = F)
write.table(MGRB_CH_hits_nofilt, "./CH_count_del_size/data/CH_hits/MGRB_CH_hits_nofilt.tsv", 
            sep="\t", row.names=F, quote=F, col.names=T)

## Plot by SEX
MGRB_CH_hits_nofilt.F <- Get_CH_hit_By_Individual_ExonicSize(MGRB_proc_SNVs_F,
                                                             MGRB_ExonicSizes.F,
                                                             filter_SNV = F)
MGRB_CH_hits_nofilt.M <- Get_CH_hit_By_Individual_ExonicSize(MGRB_proc_SNVs_M,
                                                             MGRB_ExonicSizes.M,
                                                             filter_SNV = F)
MGRB_CH_hits_nofilt.F_less20kb <- MGRB_CH_hits_nofilt.F[(which(MGRB_CH_hits_nofilt.F$exSize <= 20000)),] # restrict del size to < 20 kb
MGRB_CH_hits_nofilt.M_less20kb <- MGRB_CH_hits_nofilt.M[(which(MGRB_CH_hits_nofilt.M$exSize <= 20000)),] # restrict del size to < 20 kb

Get_CH_DelSize_Sex_Plots(MGRB_CH_hits_nofilt, 
                         CH_hits_df.F=MGRB_CH_hits_nofilt.F, 
                         CH_hits_df.M=MGRB_CH_hits_nofilt.M, "MGRB_sex")

MGRB_CH_hits_nofilt_less20kb <- MGRB_CH_hits_nofilt[(which(MGRB_CH_hits_nofilt$exSize <= 20000)),] # restrict del size to < 20 kb
Get_CH_DelSize_Sex_Plots(MGRB_CH_hits_nofilt_less20kb, 
                         CH_hits_df.F=MGRB_CH_hits_nofilt.F_less20kb, 
                         CH_hits_df.M=MGRB_CH_hits_nofilt.M_less20kb, "MGRB_sex_<20kb")


## Plot by variant type (effect_priority) - MGRB
# # MGRB nonsynonymous
# nonsyn.MGRB.SNVs <- MGRB_proc_SNVs[which(MGRB_proc_SNVs$effect_priority == "nonsynonymous SNV"),]
# nonsyn.MGRB_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(nonsyn.MGRB.SNVs, MGRB_ExonicSizes, 
#                                                             MGRB_proband_IDs, MGRB_father_IDs, MGRB_mother_IDs, 
#                                                             filter_SNV = F)
# Get_CH_DelSize_Plots(nonsyn.MGRB_CH_hits, "MGRB_nonsyn")
# 
# nonsyn.MGRB_CH_hit_less75kb <- nonsyn.MGRB_CH_hits[(which(nonsyn.MGRB_CH_hits$exSize <= 75000)),] # restrict del size to < 75 kb
# Get_CH_DelSize_Plots(nonsyn.MGRB_CH_hit_less75kb, "MGRB_nonsyn_<75kb")
# 
# # MGRB synonymous
# syn.MGRB.SNVs <- MGRB_proc_SNVs[which(MGRB_proc_SNVs$effect_priority == "synonymous SNV"),]
# syn.MGRB_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(syn.MGRB.SNVs, MGRB_ExonicSizes, 
#                                                             MGRB_proband_IDs, MGRB_father_IDs, MGRB_mother_IDs, 
#                                                             filter_SNV = F)
# Get_CH_DelSize_Plots(syn.MGRB_CH_hits, "MGRB_syn")
# 
# syn.MGRB_CH_hit_less75kb <- syn.MGRB_CH_hits[(which(syn.MGRB_CH_hits$exSize <= 75000)),] # restrict del size to < 75 kb
# Get_CH_DelSize_Plots(syn.MGRB_CH_hit_less75kb, "MGRB_syn_<75kb")



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

# #### Fisher's Test for MGRB (all variants) #### 
# MGRB_all_fisher <- Get_Fisher_Res(MGRB_CH_hits_nofilt, MGRB_CH_hits_nofilt)
# write.table(MGRB_all_fisher, "./CH_count_del_size/data/MGRB_all_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
# 
# ## MGRB+ nonsyn
# MGRB_SSC_nonsyn_fisher <- Get_Fisher_Res(MGRB_SSC_CH_hits_nofilt, nonsyn.MSSNG_SSC_CH_hits)
# write.table(MSSNG_SSC_nonsyn_fisher, "./CH_count_del_size/data/MSSNG_SSC_nonsyn_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
# 
# 
# ## MGRB syn
# MGRB_SSC_syn_fisher <- Get_Fisher_Res(MGRB_SSC_CH_hits_nofilt, syn.MGRB_SSC_CH_hits)
# write.table(MGRB_SSC_syn_fisher, "./CH_count_del_size/data/MGRB_SSC_syn_fisher.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
# 
# 

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

## MGRB
MGRB_all_fisher_sizes <- Get_Fisher_Res_SizeBins(MGRB_CH_hits_nofilt, "MGRB.all.nofilt")
MGRB_nonsyn_fisher_sizes <- Get_Fisher_Res_SizeBins(nonsyn.MGRB_CH_hits, "MGRB.nonsyn.nofilt")
MGRB_syn_fisher_sizes <- Get_Fisher_Res_SizeBins(syn.MGRB_CH_hits, "MGRB.syn.nofilt")

