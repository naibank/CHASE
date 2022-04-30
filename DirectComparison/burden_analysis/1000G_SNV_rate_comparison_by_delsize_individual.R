library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

CNV_SNVsExonicSizes_1000G <- yaml::yaml.load_file("./1000G/CNV_SNVsExonicSizes_1000G.yaml")
CNV_SNVs_1000G <- data.frame(CNV_SNVsExonicSizes_1000G[[1]])
CNV_SNVs_1000G <- CNV_SNVs_1000G[CNV_SNVs_1000G$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                       | CNV_SNVs_1000G$LoF,] 
CNV_SNVs_1000G$UID <- paste(CNV_SNVs_1000G$X.Sample, CNV_SNVs_1000G$X.id, sep='.')
CNV_SNVs_1000G <- CNV_SNVs_1000G[!(CNV_SNVs_1000G$LoF & CNV_SNVs_1000G$freq_max > 0.01), ]
CNV_ExonicSizes_1000G <- data.frame(CNV_SNVsExonicSizes_1000G[[2]])

SV_SNVsExonicSizes_1000G <- yaml::yaml.load_file("./1000G/SV_SNVsExonicSizes_1000G.yaml")
SV_SNVs_1000G <- data.frame(SV_SNVsExonicSizes_1000G[[1]])
SV_SNVs_1000G <- SV_SNVs_1000G[SV_SNVs_1000G$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                             | SV_SNVs_1000G$LoF,] 
SV_SNVs_1000G <- SV_SNVs_1000G[SV_SNVs_1000G$cnvENDAnn - SV_SNVs_1000G$cnvSTARTAnn <= 10000, ]
SV_SNVs_1000G$UID <- paste(SV_SNVs_1000G$X.Sample, SV_SNVs_1000G$X.id, sep='.')
SV_SNVs_1000G <- SV_SNVs_1000G[!(SV_SNVs_1000G$LoF & SV_SNVs_1000G$freq_max > 0.01), ]
SV_ExonicSizes_1000G <- data.frame(SV_SNVsExonicSizes_1000G[[2]])

#Combine
all_ExonicSizes_1000G <- merge(CNV_ExonicSizes_1000G, SV_ExonicSizes_1000G, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
all_ExonicSizes_1000G$exonicSize.x[is.na(all_ExonicSizes_1000G$exonicSize.x)] <- 0
all_ExonicSizes_1000G$exonicSize.y[is.na(all_ExonicSizes_1000G$exonicSize.y)] <- 0
all_ExonicSizes_1000G$exonicSize <- pmax(all_ExonicSizes_1000G$exonicSize.x, all_ExonicSizes_1000G$exonicSize.y)
combined_SNVs_1000G <- rbind(CNV_SNVs_1000G, SV_SNVs_1000G)
proc_SNVs_1000G <- dplyr::distinct(combined_SNVs_1000G, UID, .keep_all=T)

# Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
SVs_1000G <- data.frame(SV_SNVsExonicSizes_1000G[[3]])
SVs_1000G <- SVs_1000G[SVs_1000G$length <= 10000, ]
homozyg_SVs_1000G <- SVs_1000G[SVs_1000G$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
homozyg_SVs_1000G$SVUID <- with(homozyg_SVs_1000G, paste0(sample,CHROM,START,END))
proc_SNVs_1000G <- proc_SNVs_1000G[!with(proc_SNVs_1000G,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% homozyg_SVs_1000G$SVUID,]

proc_SNVs_1000G$Min_Dist <- by(proc_SNVs_1000G, seq_len(nrow(proc_SNVs_1000G)),
                                              function(r) r$MIN = min(abs(r$POS - proc_SNVs_1000G$POS[proc_SNVs_1000G$X.Sample == r$X.Sample & 
                                                                                                        proc_SNVs_1000G$CHROM == r$CHROM &
                                                                                                        proc_SNVs_1000G$POS != r$POS])) > 50)
proc_SNVs_1000G <- proc_SNVs_1000G[proc_SNVs_1000G$Min_Dist,]

Get_CH_hit_By_Individual_ExonicSize <- function(SNVs,
                                                exonic_sizes, 
                                                cumulative = F,
                                                filter_SNV = F) {
  
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

CH_hits_1000G <- Get_CH_hit_By_Individual_ExonicSize(proc_SNVs_1000G,
                                                     all_ExonicSizes_1000G)

ggplot(data=CH_hits_1000G[CH_hits_1000G$exSize <= 10000,], aes(x=exSize, y=CH_hit, group=1)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_abline(slope=sum(CH_hits_1000G$CH_hit)/sum(CH_hits_1000G$exSize), intercept=0,
              linetype='dashed',
              color='red') +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/1000G_CH_hits_smoothed_10k.png')

ggplot(data=CH_hits_1000G[CH_hits_1000G$exSize <= 10000,], aes(x=exSize, y=CH_hit)) +
  geom_smooth(aes(colour='GAM'), se=F) +
  geom_smooth(method='lm', aes(colour='Linear'), se=F) +
  geom_smooth(method='loess', aes(colour='Local Reg (loess)'), se=F) +
  geom_smooth(method='lm',formula=y~poly(x,4), aes(colour='Polynomial Order 4'), se=F) +
  geom_point() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  scale_colour_manual(name='legend', values=c('blue','red', 'brown','green')) + 
  labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
ggsave('../DT/SNV Rate Line Plots/Individual/1000G_CH_hits_smoothcomparison_10k.png')

