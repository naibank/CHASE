library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

#Load("../Unbiased/unbiased_synchevent_covariate.RData")

# Find SNV rate per every 50bp interval and plot
Get_CH_hit_By_Individual_ExonicSize <- function(SNVs,
                                                exonic_sizes, 
                                                cumulative = F,
                                                filter_SNV = T) {
  
  if (filter_SNV) {
    SNVs <- SNVs[SNVs$gnomAD_pRec >= 0.9 & SNVs$gnomAD_oe_lof_upper >= 0.35,]
  }

  df <- exonic_sizes %>% group_by(Sample.ID)
  if (!cumulative) {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% Sample.ID,]),exSize=sum(exonicSize))
  }
  else {
    df <- df %>% summarise(CH_hit=nrow(SNVs[SNVs$X.Sample %in% df$Sample.ID[df$cuts <= cuts[[1]]],])/sum(df$exonicSize[df$cuts <= cuts[[1]]]))
  }
  
  df <- data.frame(df)
  df <- df[!is.nan(df$CH_hit),]
  return (df)
}

MSSNG_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs,
                                                         MSSNG_parent_proband_all_ExonicSizes)
# MSSNG_CH_hits <- MSSNG_CH_hits[MSSNG_CH_hits$CH_hit > 0,]
ggplot(data=MSSNG_CH_hits, aes(x=exSize, y=CH_hit, group=1)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_CH_hits.png')

ggplot(data=MSSNG_CH_hits, aes(x=exSize, y=CH_hit, group=1)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope=sum(MSSNG_CH_hits$CH_hit)/sum(MSSNG_CH_hits$exSize)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_CH_hits_smoothed.png')

SSC_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(ssc_all_SNVs_combined_unique, 
                                                     ssc_all_ExonicSizes_combined_unique)
# SSC_CH_hits <- SSC_CH_hits[SSC_CH_hits$CH_hit > 0,]
ggplot(data=SSC_CH_hits, aes(x=exSize, y=CH_hit, group=1)) +
  geom_line() +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_CH_hits.png')

ggplot(data=SSC_CH_hits, aes(x=exSize, y=CH_hit, group=1)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope=sum(SSC_CH_hits$CH_hit)/sum(SSC_CH_hits$exSize)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_CH_hits_smoothed.png')

#### CNVs vs SVs in 1k to 10k range ####
MSSNG_SNV_CNV_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                           MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_MSSNG$X.Sample,])
MSSNG_SNV_CNV_CH_hits$type <- 'CNV'
# MSSNG_SNV_CNV_rates <- MSSNG_SNV_CNV_rates[MSSNG_SNV_CNV_rates$CH_hit > 0,]
MSSNG_SNV_SV_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                          MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_MSSNG$X.Sample,])
MSSNG_SNV_SV_CH_hits$type <- 'SV'
# MSSNG_SNV_SV_rates <- MSSNG_SNV_SV_rates[MSSNG_SNV_SV_rates$CH_hit > 0,]

MSSNG_SNV_CNVSV_CH_hits <- rbind(MSSNG_SNV_CNV_CH_hits,MSSNG_SNV_SV_CH_hits)
ggplot(data=MSSNG_SNV_CNVSV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_Caller_CH_hits.png')

ggplot(data=MSSNG_SNV_CNVSV_CH_hits[1000 <= MSSNG_SNV_CNVSV_CH_hits$exSize 
                                    & MSSNG_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit, group=type)) +
  geom_point(aes(color=type)) + 
  geom_smooth(method='lm',aes(color=type)) +
  geom_abline(aes(slope=sum(MSSNG_SNV_CNV_CH_hits$CH_hit)/sum(MSSNG_SNV_CNV_CH_hits$exSize),
                  intercept=0,
                  colour='CNV'),
              linetype='dashed') +
  geom_abline(aes(slope=sum(MSSNG_SNV_SV_CH_hits$CH_hit)/sum(MSSNG_SNV_SV_CH_hits$exSize),
                  intercept=0,
                  colour='SV'),
              linetype='dashed') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_Caller_CH_hits_smoothed.png')

SSC_SNV_CNV_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                         ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample,])
SSC_SNV_CNV_CH_hits$type <- 'CNV'
# SSC_SNV_CNV_CH_hits <- SSC_SNV_CNV_CH_hits[SSC_SNV_CNV_CH_hits$CH_hit > 0,]
SSC_SNV_SV_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_SV_only_SSC$X.Sample,],
                                                        ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample,])
SSC_SNV_SV_CH_hits$type <- 'SV'
# SSC_SNV_SV_CH_hits <- SSC_SNV_SV_CH_hits[SSC_SNV_SV_CH_hits$CH_hit > 0,]

SSC_SNV_CNVSV_CH_hits <- rbind(SSC_SNV_CNV_CH_hits,SSC_SNV_SV_CH_hits)
ggplot(data=SSC_SNV_CNVSV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_Caller_CH_hits.png')

ggplot(data=SSC_SNV_CNVSV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_smooth(method='lm',aes(color=type)) +
  geom_abline(aes(slope=sum(SSC_SNV_CNV_CH_hits$CH_hit)/sum(SSC_SNV_CNV_CH_hits$exSize),
                  intercept=0,
                  colour='CNV'),
              linetype='dashed') +
  geom_abline(aes(slope=sum(SSC_SNV_SV_CH_hits$CH_hit)/sum(SSC_SNV_SV_CH_hits$exSize),
                  intercept=0,
                  colour='SV'),
              linetype='dashed') +  
  geom_point(aes(color=type)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_Caller_CH_hits_smoothed.png')

#### No filter on SNVs ####
MSSNG_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs,
                                                     MSSNG_parent_proband_all_ExonicSizes,
                                                     filter_SNV = F)
# MSSNG_CH_hits_nofilt <- MSSNG_CH_hits_nofilt[MSSNG_CH_hits_nofilt$CH_hit > 0,]
ggplot(data=MSSNG_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=1)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_CH_hits_nofilt.png')

ggplot(data=MSSNG_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=1)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope=sum(MSSNG_CH_hits_nofilt$CH_hit)/sum(MSSNG_CH_hits_nofilt$exSize)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_CH_hits_nofilt_smoothed.png')

SSC_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(ssc_all_SNVs_combined_unique, 
                                                   ssc_all_ExonicSizes_combined_unique,
                                                   filter_SNV = F)
# SSC_CH_hits_nofilt <- SSC_CH_hits_nofilt[SSC_CH_hits_nofilt$CH_hit > 0,]
ggplot(data=SSC_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=1)) +
  geom_line() +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_CH_hits_nofilt.png')

ggplot(data=SSC_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=1)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope=sum(SSC_CH_hits_nofilt$CH_hit)/sum(SSC_CH_hits_nofilt$exSize)) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_CH_hits_nofilt_smoothed.png')

# Comparison of CNVs vs SVs in 1k to 10k range separately
MSSNG_SNV_CNV_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                             MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                             filter_SNV = F)
MSSNG_SNV_CNV_CH_hits_nofilt$type <- 'CNV'
# MSSNG_SNV_CNV_rates <- MSSNG_SNV_CNV_rates[MSSNG_SNV_CNV_rates$CH_hit > 0,]
MSSNG_SNV_SV_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                            MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                            filter_SNV = F)
MSSNG_SNV_SV_CH_hits_nofilt$type <- 'SV'
# MSSNG_SNV_SV_rates <- MSSNG_SNV_SV_rates[MSSNG_SNV_SV_rates$CH_hit > 0,]

MSSNG_SNV_CNVSV_CH_hits_nofilt <- rbind(MSSNG_SNV_CNV_CH_hits_nofilt,MSSNG_SNV_SV_CH_hits_nofilt)
ggplot(data=MSSNG_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_Caller_CH_hits_nofilt.png')

ggplot(data=MSSNG_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
  geom_point(aes(color=type)) + 
  geom_smooth(method='lm',aes(color=type)) +
  geom_abline(aes(slope=sum(MSSNG_SNV_CNV_CH_hits_nofilt$CH_hit)/sum(MSSNG_SNV_CNV_CH_hits_nofilt$exSize),
                  intercept=0,
                  colour='CNV'),
              linetype='dashed') +
  geom_abline(aes(slope=sum(MSSNG_SNV_SV_CH_hits_nofilt$CH_hit)/sum(MSSNG_SNV_SV_CH_hits_nofilt$exSize),
                  intercept=0,
                  colour='SV'),
              linetype='dashed') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNG_Caller_CH_hits_nofilt_smoothed.png')

SSC_SNV_CNV_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                           ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                           filter_SNV = F)
SSC_SNV_CNV_CH_hits_nofilt$type <- 'CNV'
# SSC_SNV_CNV_CH_hits_nofilt <- SSC_SNV_CNV_CH_hits_nofilt[SSC_SNV_CNV_CH_hits_nofilt$CH_hit > 0,]
SSC_SNV_SV_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_SV_only_SSC$X.Sample,],
                                                          ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample,],
                                                          filter_SNV = F)
SSC_SNV_SV_CH_hits_nofilt$type <- 'SV'
# SSC_SNV_SV_CH_hits_nofilt <- SSC_SNV_SV_CH_hits_nofilt[SSC_SNV_SV_CH_hits_nofilt$CH_hit > 0,]

SSC_SNV_CNVSV_CH_hits_nofilt <- rbind(SSC_SNV_CNV_CH_hits_nofilt,SSC_SNV_SV_CH_hits_nofilt)
ggplot(data=SSC_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_Caller_CH_hits_nofilt.png')

ggplot(data=SSC_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
  geom_smooth(method='lm',aes(color=type)) +
  geom_abline(aes(slope=sum(SSC_SNV_CNV_CH_hits_nofilt$CH_hit)/sum(SSC_SNV_CNV_CH_hits_nofilt$exSize),
                  intercept=0,
                  colour='CNV'),
              linetype='dashed') +
  geom_abline(aes(slope=sum(SSC_SNV_SV_CH_hits_nofilt$CH_hit)/sum(SSC_SNV_SV_CH_hits_nofilt$exSize),
                  intercept=0,
                  colour='SV'),
              linetype='dashed') +  
  geom_point(aes(color=type)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/SSC_Caller_CH_hits_nofilt_smoothed.png')

#### Combined MSSNG and SSC datasets ####
MSSNGSSC_SNV_CNV_CH_hits <- rbind(MSSNG_SNV_CNV_CH_hits, SSC_SNV_CNV_CH_hits)
ggplot(data=MSSNGSSC_SNV_CNV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_smooth(method='lm') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
                  intercept=0),
              linetype='dashed') +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CNVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_SV_CH_hits <- rbind(MSSNG_SNV_SV_CH_hits, SSC_SNV_SV_CH_hits)
ggplot(data=MSSNGSSC_SNV_SV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_smooth(method='lm') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
                  intercept=0),
              linetype='dashed') +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_SVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_CNVSV_CH_hits <- rbind(MSSNGSSC_SNV_CNV_CH_hits,MSSNGSSC_SNV_SV_CH_hits)
ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits[1000 <= MSSNGSSC_SNV_CNVSV_CH_hits$exSize 
                                    & MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit, group=type)) +
  geom_point(aes(color=type)) + 
  geom_smooth(method='lm',aes(color=type)) +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
                  intercept=0,
                  colour='CNV'),
              linetype='dashed') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
                  intercept=0,
                  colour='SV'),
              linetype='dashed') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CNVvsSVCaller_CH_hits_smoothed.png')


ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit)) +
  geom_point() + 
  geom_smooth(method='lm') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize),
                  intercept=0),
              linetype='dashed') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CNVSV_CH_hits_smoothed_10k.png')

#### TODO: Combined MSSNG and SSC datasets no filter####
MSSNGSSC_SNV_CNV_CH_hits <- rbind(MSSNG_SNV_CNV_CH_hits, SSC_SNV_CNV_CH_hits)
ggplot(data=MSSNGSSC_SNV_CNV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_smooth(method='lm') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
                  intercept=0),
              linetype='dashed') +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CNVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_SV_CH_hits <- rbind(MSSNG_SNV_SV_CH_hits, SSC_SNV_SV_CH_hits)
ggplot(data=MSSNGSSC_SNV_SV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
  geom_smooth(method='lm') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
                  intercept=0),
              linetype='dashed') +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_SVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_CNVSV_CH_hits <- rbind(MSSNGSSC_SNV_CNV_CH_hits,MSSNGSSC_SNV_SV_CH_hits)
ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits[1000 <= MSSNGSSC_SNV_CNVSV_CH_hits$exSize 
                                       & MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit, group=type)) +
  geom_point(aes(color=type)) + 
  geom_smooth(method='lm',aes(color=type)) +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
                  intercept=0,
                  colour='CNV'),
              linetype='dashed') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
                  intercept=0,
                  colour='SV'),
              linetype='dashed') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CNVvsSVCaller_CH_hits_smoothed.png')

#### Other nonlinear fits ####
ggplot(data=MSSNGSSC_SNV_CNV_CH_hits, aes(x=exSize, y=CH_hit)) +
  geom_smooth(aes(colour='GAM'), se=F) +
  geom_smooth(method='lm', aes(colour='Linear'), se=F) +
  geom_smooth(method='loess', aes(colour='Local Reg (loess)'), se=F) +
  geom_smooth(method='lm',formula=y~poly(x,4), aes(colour='Polynomial Order 4'), se=F) +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  scale_colour_manual(name='legend', values=c('blue','red', 'brown','green')) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CNVCaller_CH_hits_smoothcomparison.png')
#### Delineate by ASD status ####
MSSNGSSC_SNV_CNVSV_CH_hits$Relation <- NA
MSSNGSSC_SNV_CNVSV_CH_hits$Relation[MSSNGSSC_SNV_CNVSV_CH_hits$Sample.ID %in% c(MSSNG_proband_IDs,SSC_proband_IDs)] <- 'Proband'
MSSNGSSC_SNV_CNVSV_CH_hits$Relation[MSSNGSSC_SNV_CNVSV_CH_hits$Sample.ID %in% c(MSSNG_father_IDs,SSC_father_IDs)] <- 'Father'
MSSNGSSC_SNV_CNVSV_CH_hits$Relation[MSSNGSSC_SNV_CNVSV_CH_hits$Sample.ID %in% c(MSSNG_mother_IDs,SSC_mother_IDs)] <- 'Mother'
MSSNGSSC_SNV_CNVSV_CH_hits <- MSSNGSSC_SNV_CNVSV_CH_hits[!is.na(MSSNGSSC_SNV_CNVSV_CH_hits$Relation),]
ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits, aes(x=exSize, y=CH_hit)) +
  geom_smooth(method='lm') +
  geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize),
                  intercept=0),
              linetype='dashed') +
  geom_point(aes(color=Relation)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/MSSNGSSC_CH_hits_smoothed_ASD.png')

MSSNGSSC_SNV_CH_hits_probands <- MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$Relation == 'Proband',]
proband_resi <- MSSNGSSC_SNV_CH_hits_probands$CH_hit - ((sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize)) * MSSNGSSC_SNV_CH_hits_probands$exSize)
proband_resi <- proband_resi[!is.na(proband_resi)]
proband_above_expect <- sum(proband_resi > 0)
proband_below_expect <- sum(proband_resi < 0)

MSSNGSSC_SNV_CH_hits_controls <- MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$Relation %in% c('Father','Mother'),]
control_resi <- MSSNGSSC_SNV_CH_hits_controls$CH_hit - ((sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize)) * MSSNGSSC_SNV_CH_hits_controls$exSize)
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

MSSNGSSC_SNV_CH_hits_probands_less10k <- MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$Relation == 'Proband' & MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,]
proband_less10k_resi <- MSSNGSSC_SNV_CH_hits_probands_less10k$CH_hit - ((sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize)) * MSSNGSSC_SNV_CH_hits_probands_less10k$exSize)
proband_less10k_resi <- proband_less10k_resi[!is.na(proband_less10k_resi)]
proband_less10k_above_expect <- sum(proband_less10k_resi > 0)
proband_less10k_below_expect <- sum(proband_less10k_resi < 0)

MSSNGSSC_SNV_CH_hits_controls_less10k <- MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$Relation %in% c('Father','Mother') & MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,]
control_less10k_resi <- MSSNGSSC_SNV_CH_hits_controls$CH_hit - ((sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize)) * MSSNGSSC_SNV_CH_hits_controls_less10k$exSize)
control_less10k_resi <- control_less10k_resi[!is.na(control_less10k_resi)]
control_less10k_above_expect <- sum(control_less10k_resi > 0)
control_less10k_below_expect <- sum(control_less10k_resi < 0)

ratio_case_control_less10k_above_expect <- proband_less10k_above_expect/control_less10k_above_expect
ratio_case_control_less10k_below_expect <- proband_less10k_below_expect/control_less10k_below_expect
case_control_less10k_expect_OR <- ratio_case_control_less10k_above_expect/ratio_case_control_less10k_below_expect
case_control_less10k_fisher_df <- data.frame(Above=c(proband_less10k_above_expect,control_less10k_above_expect),
                                     Below=c(proband_less10k_below_expect,control_less10k_below_expect))
rownames(case_control_less10k_fisher_df) <- c('Case','Control')
case_control_less10k_fisher_res <- fisher.test(case_control_less10k_fisher_df, alternative='greater')


MSSNGSSC_SNV_CH_hits_probands_greater10k <- MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$Relation == 'Proband' & MSSNGSSC_SNV_CNVSV_CH_hits$exSize > 10000,]
proband_greater10k_resi <- MSSNGSSC_SNV_CH_hits_probands_greater10k$CH_hit - ((sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize)) * MSSNGSSC_SNV_CH_hits_probands_greater10k$exSize)
proband_greater10k_resi <- proband_greater10k_resi[!is.na(proband_greater10k_resi)]
proband_greater10k_above_expect <- sum(proband_greater10k_resi > 0)
proband_greater10k_below_expect <- sum(proband_greater10k_resi < 0)

MSSNGSSC_SNV_CH_hits_controls_greater10k <- MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$Relation %in% c('Father','Mother') & MSSNGSSC_SNV_CNVSV_CH_hits$exSize > 10000,]
control_greater10k_resi <- MSSNGSSC_SNV_CH_hits_controls$CH_hit - ((sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize)) * MSSNGSSC_SNV_CH_hits_controls_greater10k$exSize)
control_greater10k_resi <- control_greater10k_resi[!is.na(control_greater10k_resi)]
control_greater10k_above_expect <- sum(control_greater10k_resi > 0)
control_greater10k_below_expect <- sum(control_greater10k_resi < 0)

ratio_case_control_greater10k_above_expect <- proband_greater10k_above_expect/control_greater10k_above_expect
ratio_case_control_greater10k_below_expect <- proband_greater10k_below_expect/control_greater10k_below_expect
case_control_greater10k_expect_OR <- ratio_case_control_greater10k_above_expect/ratio_case_control_greater10k_below_expect
case_control_greater10k_fisher_df <- data.frame(Above=c(proband_greater10k_above_expect,control_greater10k_above_expect),
                                             Below=c(proband_greater10k_below_expect,control_greater10k_below_expect))
rownames(case_control_greater10k_fisher_df) <- c('Case','Control')
case_control_greater10k_fisher_res <- fisher.test(case_control_greater10k_fisher_df, alternative='greater')



