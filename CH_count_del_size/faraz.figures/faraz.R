
#SNV Counts for all CNV
SNV_in_CNV_only_MSSNG <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$UID %in% MSSNG_parent_proband_SNVs$UID, ]
SNV_in_CNV_only_SSC <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$UID %in% SSC_parent_proband_SNVs$UID,]



MSSNG_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(MSSNG_parent_proband_proc_SNVs,
                                                            MSSNG_parent_proband_all_ExonicSizes,
                                                            filter_SNV = F)
MSSNG_CH_hits_nofilt <- MSSNG_CH_hits_nofilt[MSSNG_CH_hits_nofilt$CH_hit > 0,]
# 
# ggplot(data=MSSNG_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=1)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   geom_abline(slope=sum(MSSNG_CH_hits_nofilt$CH_hit)/sum(MSSNG_CH_hits_nofilt$exSize)) +
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNG_CH_hits_nofilt_smoothed.png')

SSC_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs, 
                                                          SSC_parent_proband_all_ExonicSizes,
                                                          filter_SNV = F)
SSC_CH_hits_nofilt <- SSC_CH_hits_nofilt[SSC_CH_hits_nofilt$CH_hit > 0,]

# ggplot(data=SSC_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=1)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   geom_abline(slope=sum(SSC_CH_hits_nofilt$CH_hit)/sum(SSC_CH_hits_nofilt$exSize)) +
#   theme(axis.text.x=element_text(angle=90, hjust=1)) +
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/SSC_CH_hits_nofilt_smoothed.png')

#####################################################################################################################################################################################################################################################

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



SSC_SNV_CNV_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$X.Sample %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                           SSC_parent_proband_all_ExonicSizes[SSC_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample,])
SSC_SNV_CNV_CH_hits$type <- 'CNV'
# SSC_SNV_CNV_CH_hits <- SSC_SNV_CNV_CH_hits[SSC_SNV_CNV_CH_hits$CH_hit > 0,]
SSC_SNV_SV_CH_hits <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$X.Sample %in% SNV_in_SV_only_SSC$X.Sample,],
                                                          SSC_parent_proband_all_ExonicSizes[SSC_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample,])
SSC_SNV_SV_CH_hits$type <- 'SV'
# SSC_SNV_SV_CH_hits <- SSC_SNV_SV_CH_hits[SSC_SNV_SV_CH_hits$CH_hit > 0,]

SSC_SNV_CNVSV_CH_hits <- rbind(SSC_SNV_CNV_CH_hits,SSC_SNV_SV_CH_hits)


MSSNG_SNV_CNVSV_CH_hits_nofilt <- rbind(MSSNG_SNV_CNV_CH_hits_nofilt,MSSNG_SNV_SV_CH_hits_nofilt)

# ggplot(data=MSSNG_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_line(aes(color=type)) +
#   geom_point(aes(color=type)) + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNG_Caller_CH_hits_nofilt.png')
# 
# ggplot(data=MSSNG_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_point(aes(color=type)) + 
#   geom_smooth(method='lm',aes(color=type)) +
#   geom_abline(aes(slope=sum(MSSNG_SNV_CNV_CH_hits_nofilt$CH_hit)/sum(MSSNG_SNV_CNV_CH_hits_nofilt$exSize),
#                   intercept=0,
#                   colour='CNV'),
#               linetype='dashed') +
#   geom_abline(aes(slope=sum(MSSNG_SNV_SV_CH_hits_nofilt$CH_hit)/sum(MSSNG_SNV_SV_CH_hits_nofilt$exSize),
#                   intercept=0,
#                   colour='SV'),
#               linetype='dashed') +
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNG_Caller_CH_hits_nofilt_smoothed.png')

SSC_SNV_CNV_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$X.Sample %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                                  SSC_parent_proband_all_ExonicSizes[SSC_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample,],
                                                                  filter_SNV = F)
SSC_SNV_CNV_CH_hits_nofilt$type <- 'CNV'
# SSC_SNV_CNV_CH_hits_nofilt <- SSC_SNV_CNV_CH_hits_nofilt[SSC_SNV_CNV_CH_hits_nofilt$CH_hit > 0,]
SSC_SNV_SV_CH_hits_nofilt <- Get_CH_hit_By_Individual_ExonicSize(SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$X.Sample %in% SNV_in_SV_only_SSC$X.Sample,],
                                                                 SSC_parent_proband_all_ExonicSizes[SSC_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample,],
                                                                 filter_SNV = F)
SSC_SNV_SV_CH_hits_nofilt$type <- 'SV'
# SSC_SNV_SV_CH_hits_nofilt <- SSC_SNV_SV_CH_hits_nofilt[SSC_SNV_SV_CH_hits_nofilt$CH_hit > 0,]

SSC_SNV_CNVSV_CH_hits_nofilt <- rbind(SSC_SNV_CNV_CH_hits_nofilt,SSC_SNV_SV_CH_hits_nofilt)
# ggplot(data=SSC_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_line(aes(color=type)) +
#   geom_point(aes(color=type)) + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/SSC_Caller_CH_hits_nofilt.png')
# 
# ggplot(data=SSC_SNV_CNVSV_CH_hits_nofilt, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_smooth(method='lm',aes(color=type)) +
#   geom_abline(aes(slope=sum(SSC_SNV_CNV_CH_hits_nofilt$CH_hit)/sum(SSC_SNV_CNV_CH_hits_nofilt$exSize),
#                   intercept=0,
#                   colour='CNV'),
#               linetype='dashed') +
#   geom_abline(aes(slope=sum(SSC_SNV_SV_CH_hits_nofilt$CH_hit)/sum(SSC_SNV_SV_CH_hits_nofilt$exSize),
#                   intercept=0,
#                   colour='SV'),
#               linetype='dashed') +  
#   geom_point(aes(color=type)) + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/SSC_Caller_CH_hits_nofilt_smoothed.png')

#### Combined MSSNG and SSC datasets ####
MSSNGSSC_SNV_CNV_CH_hits <- rbind(MSSNG_SNV_CNV_CH_hits, SSC_SNV_CNV_CH_hits)
# ggplot(data=MSSNGSSC_SNV_CNV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_smooth(method='lm') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
#                   intercept=0),
#               linetype='dashed') +
#   geom_point() + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_CNVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_SV_CH_hits <- rbind(MSSNG_SNV_SV_CH_hits, SSC_SNV_SV_CH_hits)
# ggplot(data=MSSNGSSC_SNV_SV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_smooth(method='lm') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
#                   intercept=0),
#               linetype='dashed') +
#   geom_point() + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_SVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_CNVSV_CH_hits <- rbind(MSSNGSSC_SNV_CNV_CH_hits,MSSNGSSC_SNV_SV_CH_hits)
# ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits[1000 <= MSSNGSSC_SNV_CNVSV_CH_hits$exSize 
#                                        & MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit, group=type)) +
#   geom_point(aes(color=type)) + 
#   geom_smooth(method='lm',aes(color=type)) +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
#                   intercept=0,
#                   colour='CNV'),
#               linetype='dashed') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
#                   intercept=0,
#                   colour='SV'),
#               linetype='dashed') +
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_CNVvsSVCaller_CH_hits_smoothed.png')
# 
# 
# ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits[MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit)) +
#   geom_point() + 
#   geom_smooth(method='lm') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize),
#                   intercept=0),
#               linetype='dashed') +
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_CNVSV_CH_hits_smoothed_10k.png')

#### TODO: Combined MSSNG and SSC datasets no filter####
MSSNGSSC_SNV_CNV_CH_hits <- rbind(MSSNG_SNV_CNV_CH_hits, SSC_SNV_CNV_CH_hits)
# ggplot(data=MSSNGSSC_SNV_CNV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_smooth(method='lm') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
#                   intercept=0),
#               linetype='dashed') +
#   geom_point() + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_CNVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_SV_CH_hits <- rbind(MSSNG_SNV_SV_CH_hits, SSC_SNV_SV_CH_hits)
# ggplot(data=MSSNGSSC_SNV_SV_CH_hits, aes(x=exSize, y=CH_hit, group=type)) +
#   geom_smooth(method='lm') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
#                   intercept=0),
#               linetype='dashed') +
#   geom_point() + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_SVCaller_CH_hits_smoothed.png')

MSSNGSSC_SNV_CNVSV_CH_hits <- rbind(MSSNGSSC_SNV_CNV_CH_hits,MSSNGSSC_SNV_SV_CH_hits)
# ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits[1000 <= MSSNGSSC_SNV_CNVSV_CH_hits$exSize 
#                                        & MSSNGSSC_SNV_CNVSV_CH_hits$exSize <= 10000,], aes(x=exSize, y=CH_hit, group=type)) +
#   geom_point(aes(color=type)) + 
#   geom_smooth(method='lm',aes(color=type)) +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNV_CH_hits$exSize),
#                   intercept=0,
#                   colour='CNV'),
#               linetype='dashed') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_SV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_SV_CH_hits$exSize),
#                   intercept=0,
#                   colour='SV'),
#               linetype='dashed') +
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_CNVvsSVCaller_CH_hits_smoothed.png')


#### Delineate by ASD status ####
MSSNGSSC_SNV_CNVSV_CH_hits$Relation <- NA
MSSNGSSC_SNV_CNVSV_CH_hits$Relation[MSSNGSSC_SNV_CNVSV_CH_hits$Sample.ID %in% c(MSSNG_proband_IDs,SSC_proband_IDs)] <- 'Proband'
MSSNGSSC_SNV_CNVSV_CH_hits$Relation[MSSNGSSC_SNV_CNVSV_CH_hits$Sample.ID %in% c(MSSNG_father_IDs,SSC_father_IDs)] <- 'Father'
MSSNGSSC_SNV_CNVSV_CH_hits$Relation[MSSNGSSC_SNV_CNVSV_CH_hits$Sample.ID %in% c(MSSNG_mother_IDs,SSC_mother_IDs)] <- 'Mother'
MSSNGSSC_SNV_CNVSV_CH_hits <- MSSNGSSC_SNV_CNVSV_CH_hits[!is.na(MSSNGSSC_SNV_CNVSV_CH_hits$Relation),]
# ggplot(data=MSSNGSSC_SNV_CNVSV_CH_hits, aes(x=exSize, y=CH_hit)) +
#   geom_smooth(method='lm') +
#   geom_abline(aes(slope=sum(MSSNGSSC_SNV_CNVSV_CH_hits$CH_hit)/sum(MSSNGSSC_SNV_CNVSV_CH_hits$exSize),
#                   intercept=0),
#               linetype='dashed') +
#   geom_point(aes(color=Relation)) + 
#   theme(axis.text.x=element_text(angle=90, hjust=1)) + 
#   labs(y='No. of CH Events', x='Total deleted exonic size (bp)')
# ggsave('./CH_count_del_size/faraz.figures/MSSNGSSC_CH_hits_smoothed_ASD.png', width = 7, height = 5)


# write.table(MSSNGSSC_SNV_CNVSV_CH_hits, "./CH_count_del_size/faraz.MSSNGSSC_SNV_CNVSV_CH_hits.with.filtering.tsv", sep="\t", row.names=F, quote=F, col.names=T)
