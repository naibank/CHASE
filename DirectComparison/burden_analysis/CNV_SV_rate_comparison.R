library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

#Load("../Unbiased/unbiased_synchevent_covariate.RData")

# Find SNV rate per every 50bp interval and plot
Get_SNV_Rate_By_ExonicSize_Interval <- function(SNVs, exonic_sizes, cumulative = F,
                                                bin_width=50, from=50,to=3000000) {
  bin_breaks <- seq(from=from, to=to, by=bin_width)
  df <- exonic_sizes %>% mutate(cuts=cut(exonic_sizes$exonicSize, 
                                         breaks=bin_breaks,
                                         include.lowest=T,
                                         right=F,
                                         ordered_result=T)) %>% group_by(cuts)
  df <- df[!is.na(df$cuts),]
  if (!cumulative) {
    df <- df %>% summarise(SNV_rate=nrow(SNVs[SNVs$X.Sample %in% Sample.ID,])/sum(exonicSize))
  }
  else {
    df <- df %>% summarise(SNV_rate=nrow(SNVs[SNVs$X.Sample %in% df$Sample.ID[df$cuts <= cuts[[1]]],])/sum(df$exonicSize[df$cuts <= cuts[[1]]]))
  }
  
  df <- data.frame(df)
  return (df)
}

x_labels <- c('[50,','[1000,','[10000,')
MSSNG_SNV_rates <- Get_SNV_Rate_By_ExonicSize_Interval(MSSNG_parent_proband_proc_SNVs, 
                                             MSSNG_parent_proband_all_ExonicSizes,
                                             bin_width=50)
ggplot(data=MSSNG_SNV_rates, aes(x=cuts, y=SNV_rate, group=1)) +
  geom_line() +
  geom_point() + 
  scale_x_discrete(breaks=MSSNG_SNV_rates$cuts[c(which(Reduce('|',lapply(x_labels, startsWith, x=as.character(MSSNG_SNV_rates$cuts)))),nrow(MSSNG_SNV_rates))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/MSSNG_SNV_rates.png')

SSC_SNV_rates <- Get_SNV_Rate_By_ExonicSize_Interval(ssc_all_SNVs_combined_unique, 
                                                     ssc_all_ExonicSizes_combined_unique,
                                                     bin_width=50)
ggplot(data=SSC_SNV_rates, aes(x=cuts, y=SNV_rate, group=1)) +
  geom_line() +
  geom_point() + 
  scale_x_discrete(breaks=SSC_SNV_rates$cuts[c(which(Reduce('|',lapply(x_labels, startsWith, x=as.character(SSC_SNV_rates$cuts)))),nrow(SSC_SNV_rates))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/SSC_SNV_rates.png')

# Comparison of CNVs vs SVs in 1k to 10k range separately
MSSNG_SNV_CNV_rates <- Get_SNV_Rate_By_ExonicSize_Interval(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                           MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                           bin_width=50,
                                                           from=1000,
                                                           to=10000)
MSSNG_SNV_CNV_rates$type <- 'CNV'
MSSNG_SNV_SV_rates <- Get_SNV_Rate_By_ExonicSize_Interval(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                           MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                           bin_width=50,
                                                           from=1000,
                                                           to=10000)
MSSNG_SNV_SV_rates$type <- 'SV'
ggplot(data=rbind(MSSNG_SNV_CNV_rates,MSSNG_SNV_SV_rates), aes(x=cuts, y=SNV_rate, group=type)) +
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) + 
  scale_x_discrete(breaks=rbind(MSSNG_SNV_CNV_rates,MSSNG_SNV_SV_rates)$cuts[c(1,nrow(MSSNG_SNV_CNV_rates))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/MSSNG_1kto10k_SNV_rates.png')
SSC_SNV_CNV_rates <- Get_SNV_Rate_By_ExonicSize_Interval(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                         ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                         bin_width=50,
                                                         from=1000,
                                                         to=10000)
SSC_SNV_CNV_rates$type <- 'CNV'
SSC_SNV_SV_rates <- Get_SNV_Rate_By_ExonicSize_Interval(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_SV_only_SSC$X.Sample,],
                                                         ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample,],
                                                         bin_width=50,
                                                         from=1000,
                                                         to=10000)
SSC_SNV_SV_rates$type <- 'SV'
ggplot(data=rbind(SSC_SNV_CNV_rates,SSC_SNV_SV_rates), aes(x=cuts, y=SNV_rate, group=type)) +
  geom_line(aes(color=type)) +
  geom_point(aes(color=type)) + 
  scale_x_discrete(breaks=rbind(SSC_SNV_CNV_rates,SSC_SNV_SV_rates)$cuts[c(1,nrow(SSC_SNV_CNV_rates))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/SSC_1kto10k_SNV_rates.png')


#### Cumulative ####
MSSNG_SNV_rates_cum <- Get_SNV_Rate_By_ExonicSize_Interval(MSSNG_parent_proband_proc_SNVs, 
                                                       MSSNG_parent_proband_all_ExonicSizes,
                                                       bin_width=50,
                                                       cumulative=T)
ggplot(data=MSSNG_SNV_rates_cum, aes(x=cuts, y=SNV_rate, group=1)) +
  geom_line() +
  scale_x_discrete(breaks=MSSNG_SNV_rates_cum$cuts[c(which(Reduce('|',lapply(x_labels, startsWith, x=as.character(MSSNG_SNV_rates_cum$cuts)))),nrow(MSSNG_SNV_rates_cum))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/MSSNG_SNV_rates_cum.png')

SSC_SNV_rates_cum <- Get_SNV_Rate_By_ExonicSize_Interval(ssc_all_SNVs_combined_unique, 
                                                     ssc_all_ExonicSizes_combined_unique,
                                                     bin_width=50,
                                                     cumulative=T)
ggplot(data=SSC_SNV_rates_cum[-1:-5,], aes(x=cuts, y=SNV_rate, group=1)) +
  geom_line() +
  scale_x_discrete(breaks=SSC_SNV_rates_cum$cuts[c(which(Reduce('|',lapply(x_labels, startsWith, x=as.character(SSC_SNV_rates_cum$cuts)))),nrow(SSC_SNV_rates_cum))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/SSC_SNV_rates_cum_f5r.png')

# Comparison of CNVs vs SVs in 1k to 10k range separately
MSSNG_SNV_CNV_rates_cum <- Get_SNV_Rate_By_ExonicSize_Interval(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                           MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_MSSNG$X.Sample,],
                                                           bin_width=50,
                                                           from=1000,
                                                           to=10000,
                                                           cumulative=T)
MSSNG_SNV_CNV_rates_cum$type <- 'CNV'
MSSNG_SNV_SV_rates_cum <- Get_SNV_Rate_By_ExonicSize_Interval(MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                          MSSNG_parent_proband_all_ExonicSizes[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_MSSNG$X.Sample,],
                                                          bin_width=50,
                                                          from=1000,
                                                          to=10000,
                                                          cumulative=T)
MSSNG_SNV_SV_rates_cum$type <- 'SV'
ggplot(data=rbind(MSSNG_SNV_CNV_rates_cum,MSSNG_SNV_SV_rates_cum), aes(x=cuts, y=SNV_rate, group=type)) +
  geom_line(aes(color=type)) +
  scale_x_discrete(breaks=rbind(MSSNG_SNV_CNV_rates_cum,MSSNG_SNV_SV_rates_cum)$cuts[c(1,nrow(MSSNG_SNV_CNV_rates_cum))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/MSSNG_1kto10k_SNV_rates_cum.png')

SSC_SNV_CNV_rates_cum <- Get_SNV_Rate_By_ExonicSize_Interval(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                         ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                         bin_width=50,
                                                         from=1000,
                                                         to=10000,
                                                         cumulative=T)
SSC_SNV_CNV_rates_cum$type <- 'CNV'
SSC_SNV_SV_rates_cum <- Get_SNV_Rate_By_ExonicSize_Interval(ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$X.Sample %in% SNV_in_SV_only_SSC$X.Sample,],
                                                        ssc_all_ExonicSizes_combined_unique[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample,],
                                                        bin_width=50,
                                                        from=1000,
                                                        to=10000,
                                                        cumulative=T)
SSC_SNV_SV_rates_cum$type <- 'SV'
ggplot(data=rbind(SSC_SNV_CNV_rates_cum,SSC_SNV_SV_rates_cum), aes(x=cuts, y=SNV_rate, group=type)) +
  geom_line(aes(color=type)) +
  scale_x_discrete(breaks=rbind(SSC_SNV_CNV_rates_cum,SSC_SNV_SV_rates_cum)$cuts[c(1,nrow(SSC_SNV_CNV_rates_cum))]) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='SNV rate per bp', x='bp Bin')
ggsave('../DT/SNV Rate Line Plots/SSC_1kto10k_SNV_rates_cum.png')

# Get SVs to check for false calls
MSSNG_false_SV_check <- MSSNG_parent_proband_SVs_ILMN
MSSNG_false_SV_check$ID <- paste(MSSNG_false_SV_check$chrAnn,MSSNG_false_SV_check$START,MSSNG_false_SV_check$END)
MSSNG_false_SV_check <- MSSNG_false_SV_check[MSSNG_false_SV_check$ID %in% paste(MSSNG_parent_proband_proc_SNVs$cnvchrAnn,
                                                              MSSNG_parent_proband_proc_SNVs$cnvSTARTAnn,
                                                              MSSNG_parent_proband_proc_SNVs$cnvENDAnn),]
MSSNG_false_SV_check <- MSSNG_false_SV_check[1000 <= MSSNG_false_SV_check$length 
                                 & MSSNG_false_SV_check$length <= 2500,]
MSSNG_false_SV_check <- MSSNG_false_SV_check %>% distinct(ID, .keep_all=T) %>% arrange(length)
MSSNG_false_SV_check <- MSSNG_false_SV_check[, c('sample','chrAnn','START','END','length')]
write.csv(MSSNG_false_SV_check, 'False SV check/MSSNG_SV_1kto2.5k.csv')


SSC_false_SV_check <- SSC_parent_child_SVs 
SSC_false_SV_check$ID <- paste(SSC_false_SV_check$chrAnn,SSC_false_SV_check$START,SSC_false_SV_check$END)
SSC_false_SV_check <- SSC_false_SV_check[SSC_false_SV_check$ID %in% paste(ssc_all_SNVs_combined_unique$cnvchrAnn,
                                                                          ssc_all_SNVs_combined_unique$cnvSTARTAnn,
                                                                          ssc_all_SNVs_combined_unique$cnvENDAnn),]
SSC_false_SV_check <- SSC_false_SV_check[1000 <= SSC_false_SV_check$length 
                                             & SSC_false_SV_check$length <= 2500,]
SSC_false_SV_check <- SSC_false_SV_check %>% distinct(ID, .keep_all=T) %>% arrange(length)
SSC_false_SV_check <- SSC_false_SV_check[, c('sample','chrAnn','START','END','length')]
write.csv(SSC_false_SV_check, 'False SV check/SSC_SV_1kto2.5k.csv')

# Check false SNV in SV spike for MSSNG and SSC
# spike index 8 and 12 for mssng
MSSNG_false_SNVinSV_check <- SNV_in_SV_only_MSSNG[SNV_in_SV_only_MSSNG$X.Sample 
                                                  %in% 
                                                    MSSNG_parent_proband_all_ExonicSizes$Sample.ID[(MSSNG_parent_proband_all_ExonicSizes$exonicSize >= 1350 & MSSNG_parent_proband_all_ExonicSizes$exonicSize <= 1400) |
                                                                                                     (MSSNG_parent_proband_all_ExonicSizes$exonicSize >= 1550 & MSSNG_parent_proband_all_ExonicSizes$exonicSize <= 1600)],]
MSSNG_false_SNVinSV_check <- MSSNG_false_SNVinSV_check[, c('X.Sample','CHROM','POS','cnvSTARTAnn','cnvENDAnn')]
write.csv(MSSNG_false_SNVinSV_check, 'False SV check/MSSNG_SNV_falsecheck.csv')

# spike index 7 and 37,38 and 101 for ssc
SSC_false_SNVinSV_check <- SNV_in_SV_only_SSC[SNV_in_SV_only_SSC$X.Sample 
                                                  %in% 
                                                    ssc_all_ExonicSizes_combined_unique$Sample.ID[(ssc_all_ExonicSizes_combined_unique$exonicSize >= 1300 & ssc_all_ExonicSizes_combined_unique$exonicSize <= 1350) |
                                                                                                    # (ssc_all_ExonicSizes_combined_unique$exonicSize >= 2800 & ssc_all_ExonicSizes_combined_unique$exonicSize <= 2850) |
                                                                                                    (ssc_all_ExonicSizes_combined_unique$exonicSize >= 6000 & ssc_all_ExonicSizes_combined_unique$exonicSize <= 6050)],]
SSC_false_SNVinSV_check <- SSC_false_SNVinSV_check[, c('X.Sample','CHROM','POS','cnvSTARTAnn','cnvENDAnn')]
write.csv(SSC_false_SNVinSV_check, 'False SV check/SSC_SNV_falsecheck.csv')
