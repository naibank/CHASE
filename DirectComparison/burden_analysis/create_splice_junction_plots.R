library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

#Load("../Unbiased/unbiased_synchevent_covariate.RData")

Generate_Splice_Junction_Distribution_Plot <- function(meta,
                                                       SNV_data,
                                                       name,
                                                       title,
                                                       proband_IDs,
                                                       father_IDs,
                                                       mother_IDs,
                                                       type='synonymous SNV',
                                                       save_path='../DT/Splice Junction Plots/') {
  sj_df <- SNV_data[SNV_data$X.Sample %in% meta$Sample.ID
                    & SNV_data$effect_priority == type
                    & SNV_data$gnomAD_pRec >= 0.9 & !is.na(SNV_data$gnomAD_pRec)
                    & SNV_data$gnomAD_oe_lof_upper >= 0.35 & !is.na(SNV_data$gnomAD_oe_lof_upper)
                    & SNV_data$distance_spliceJunction <= 15
                    , c('X.Sample', 'distance_spliceJunction')]
  sj_df$Group <- 'NA'
  sj_df$Group[sj_df$X.Sample %in% proband_IDs] <- 'Proband'
  sj_df$Group[sj_df$X.Sample %in% father_IDs] <- 'Father'
  sj_df$Group[sj_df$X.Sample %in% mother_IDs] <- 'Mother'

  ggplot(sj_df, aes(x=Group, y=distance_spliceJunction, color=Group)) +
    geom_boxplot(notch = T, outlier.shape = NA) +
    stat_summary(fun=mean) +
    scale_x_discrete(limits=c('Proband','Father','Mother')) +
    coord_cartesian(ylim=quantile(sj_df$distance_spliceJunction, c(0.1, 0.9))) +
    coord_cartesian(ylim=c(0,15)) +
    labs(y='Distance from Splice Junction', x='Group', title=title) +
    theme(plot.title=element_text(hjust=0.5))
  
  ggsave(paste(save_path, paste(name,'.png')))
}

Generate_Splice_Junction_Distribution_Plot(MSSNG_meta,
                                           MSSNG_parent_proband_proc_SNVs,
                                           'MSSNG_parent_proband_proc_SNVs_15cutoff',
                                           "MSSNG parent-proband",
                                           MSSNG_proband_IDs,
                                           MSSNG_father_IDs,
                                           MSSNG_mother_IDs)
Generate_Splice_Junction_Distribution_Plot(MSSNG_meta,
                                           MSSNG_parent_proband_0.1CH_SNVs,
                                           'MSSNG_parent_proband_0.1CH_SNVs_15cutoff',
                                           "MSSNG parent-proband 0.1% CH",
                                           MSSNG_proband_IDs,
                                           MSSNG_father_IDs,
                                           MSSNG_mother_IDs)
Generate_Splice_Junction_Distribution_Plot(SSC_meta,
                                           SSC_parent_proband_0.01CH_SNVs,
                                           'SSC_parent_proband_0.01CH_SNVs_15cutoff',
                                           "SSC parent-proband 0.01% CH",
                                           SSC_proband_IDs,
                                           SSC_father_IDs,
                                           SSC_mother_IDs)

