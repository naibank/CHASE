library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/Volumes/T5 ExFAT/Co-op/S22 Co-op SickKids (Scherer Lab)/CHASE_stor/SPARK/SPARK_TDT/")

Load_TTRes_DF <- function(path, add_cols = F) {
  res <- fread(path, data.table=F)
  res <- res[res$event_freq == 0.01 & res$pRec == 0.9,]
  res$target_UT <- res$target.parent - res$target.child
  res$target.child <- as.character(res$target.child)
  res$target_UT <- as.character(res$target_UT)
  
  return (data.frame(res))
}

Get_Interpolated_Text_TT_DF <- function(df, n=21, skip_lof=F, is_probUS=F) {
  df$P <- -log10(df$P)
  res <- df[, c('SNV_freq', 'P', 'variant_type', 'target.child', 'target_UT')]
  res$count_text <- paste('(',res$target.child,',',res$target_UT,')',sep='')
  
  df_text_all <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$variant_type=='damaging_missense'], df$P[df$variant_type=='damaging_missense'], n=n)))
  colnames(df_text_all) <- c('SNV_freq', 'P')
  df_text_all$variant_type <- 'damaging_missense'
  df_text_all$target.child <- ""
  df_text_all$target_UT <- ""
  df_text_all$count_text <- ""
  
  df_text_syn <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$variant_type=='synonymous'], df$P[df$variant_type=='synonymous'], n=n)))
  colnames(df_text_syn) <- c('SNV_freq', 'P')
  df_text_syn$variant_type <- 'synonymous'
  df_text_syn$target.child <- ""
  df_text_syn$target_UT <- ""
  df_text_syn$count_text <- ""
  
  df_text_mis <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$variant_type=='nonsynonymous'], df$P[df$variant_type=='nonsynonymous'], n=n)))
  colnames(df_text_mis) <- c('SNV_freq', 'P')
  df_text_mis$variant_type <- 'nonsynonymous'
  df_text_mis$target.child <- ""
  df_text_mis$target_UT <- ""
  df_text_mis$count_text <- ""
  
  final <- bind_rows(res, df_text_all, df_text_syn, df_text_mis)
  
  if (!skip_lof) {
    df_text_lof <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$variant_type=='LoF'], df$P[df$variant_type=='LoF'], n=n)))
    colnames(df_text_lof) <- c('SNV_freq', 'P')
    df_text_lof$variant_type <- 'LoF'
    df_text_lof$target.child <- ""
    df_text_lof$target_UT <- ""
    df_text_lof$count_text <- ""
    final <- bind_rows(final, df_text_lof)
  }
  
  return (final)
}

Create_TT_Result_Plot <- function(df, df_text, title, name, save_path='./FarazFinalPresPlots/') {
  ggplot(data=df, aes(x=SNV_freq, y=-log10(P), color=variant_type)) +
    geom_point() +
    # geom_point(data=df[df$P <= 0.05,],
    #            aes(x=SNV_freq,y=coeff,color=variant_type),
    #            shape=9,
    #            size=3) +
    geom_line(aes(group=variant_type)) +
    scale_x_discrete(limits=c('100%', '10%', '1%')) +
    geom_hline(yintercept=-log10(0.05), linetype='dashed',color='red') +
    geom_text(aes(x=0.55, y=-log10(0.05)+0.05, label='p = 0.05'), color='red', size=2) +
    geom_text_repel(data=df_text,
                    mapping=aes(y=P, x=SNV_freq, label=count_text),
                    force=10,
                    max.overlaps = Inf,
                    size=3,
                    segment.linetype=2) +
    # geom_point(data=df_text, aes(x=SNV_freq, y=P, color=variant_type)) +
    labs(y='-log10(P)', x='SNV MAF', color="Variant Type", title=title) + 
    theme(plot.title=element_text(hjust=0.5))
  ggsave(paste(save_path, paste(name,'.png',sep='')), width=6, height=6)
}


#### SPARK parent-proband Results ####
SPARK_SNV_100P_TT_df <- Load_TTRes_DF('./data/ SPARK_proband_FisherTest_eventfreq_snvfreq1.tsv')
SPARK_SNV_100P_TT_df$SNV_freq <- 1
SPARK_SNV_10P_TT_df <- Load_TTRes_DF('./data/ SPARK_proband_FisherTest_eventfreq_snvfreq0.1.tsv')
SPARK_SNV_10P_TT_df$SNV_freq <- 2
SPARK_SNV_1P_TT_df <- Load_TTRes_DF('./data/ SPARK_proband_FisherTest_eventfreq_snvfreq0.01.tsv')
SPARK_SNV_1P_TT_df$SNV_freq <- 3

SPARK_TT_df <- rbind(SPARK_SNV_100P_TT_df,SPARK_SNV_10P_TT_df,SPARK_SNV_1P_TT_df)
SPARK_text_TT_df <- Get_Interpolated_Text_TT_DF(SPARK_TT_df, 50)

Create_TT_Result_Plot(SPARK_TT_df, SPARK_text_TT_df, 'SPARK Parent-Proband', 'SPARK_TTP')

#### SPARK parent-US Results ####
SPARKUS_SNV_100P_TT_df <- Load_TTRes_DF('./data/ SPARK_unaffected_sibling_FisherTest_eventfreq_snvfreq1.tsv')
SPARKUS_SNV_100P_TT_df$SNV_freq <- 1
SPARKUS_SNV_10P_TT_df <- Load_TTRes_DF('./data/ SPARK_unaffected_sibling_FisherTest_eventfreq_snvfreq0.1.tsv')
SPARKUS_SNV_10P_TT_df$SNV_freq <- 2
SPARKUS_SNV_1P_TT_df <- Load_TTRes_DF('./data/ SPARK_unaffected_sibling_FisherTest_eventfreq_snvfreq0.01.tsv')
SPARKUS_SNV_1P_TT_df$SNV_freq <- 3

SPARKUS_TT_df <- rbind(SPARKUS_SNV_100P_TT_df,SPARKUS_SNV_10P_TT_df,SPARKUS_SNV_1P_TT_df)
SPARKUS_text_TT_df <- Get_Interpolated_Text_TT_DF(SPARKUS_TT_df, 50)

Create_TT_Result_Plot(SPARKUS_TT_df, SPARKUS_text_TT_df, 'SPARK Parent-Unaffected Sibling', 'SPARK_US_TTP')


# 
# 
# #### MSSNG parent-proband Results ####
# MSSNG_SNV_100P_TT_df <- Load_TTRes_DF('../TT/Combined Data/MSSNG/MSSNG_FisherTest_eventfreq.tsv')
# MSSNG_SNV_100P_TT_df$SNV_freq <- 1
# MSSNG_SNV_10P_TT_df <- Load_TTRes_DF('../TT/Combined Data/MSSNG/SNV0.1/MSSNG_FisherTest_eventfreq_snvfreq0.1.tsv')
# MSSNG_SNV_10P_TT_df$SNV_freq <- 2
# MSSNG_SNV_1P_TT_df <- Load_TTRes_DF('../TT/Combined Data/MSSNG/SNV0.01/MSSNG_FisherTest_eventfreq_snvfreq0.01.tsv')
# MSSNG_SNV_1P_TT_df$SNV_freq <- 3
# 
# MSSNG_TT_df <- rbind(MSSNG_SNV_100P_TT_df,MSSNG_SNV_10P_TT_df,MSSNG_SNV_1P_TT_df)
# MSSNG_text_TT_df <- Get_Interpolated_Text_TT_DF(MSSNG_TT_df, 50)
# 
# Create_TT_Result_Plot(MSSNG_TT_df, MSSNG_text_TT_df, 'MSSNG Parent-Proband', 'MSSNG_TTP')
# 
# 
# #### SSC parent-proband Results ####
# SSC_SNV_100P_TT_df <- Load_TTRes_DF('../TT/Combined Data/SSC/ SSC_proband_FisherTest_eventfreq.tsv')
# SSC_SNV_100P_TT_df$SNV_freq <- 1
# SSC_SNV_10P_TT_df <- Load_TTRes_DF('../TT/Combined Data/SSC/SNV0.1/ SSC_proband_FisherTest_eventfreq_snvfreq0.1.tsv')
# SSC_SNV_10P_TT_df$SNV_freq <- 2
# SSC_SNV_1P_TT_df <- Load_TTRes_DF('../TT/Combined Data/SSC/SNV0.01/ SSC_proband_FisherTest_eventfreq_snvfreq0.01.tsv')
# SSC_SNV_1P_TT_df$SNV_freq <- 3
# 
# SSC_TT_df <- rbind(SSC_SNV_100P_TT_df,SSC_SNV_10P_TT_df,SSC_SNV_1P_TT_df)
# SSC_text_TT_df <- Get_Interpolated_Text_TT_DF(SSC_TT_df, 50)
# 
# Create_TT_Result_Plot(SSC_TT_df, SSC_text_TT_df, 'SSC Parent-Proband', 'SSC_TTP')
# 
# #### SSC parent-US Results ####
# SSCUS_SNV_100P_TT_df <- Load_TTRes_DF('../TT/Combined Data/SSC/ SSC_unaffected_sibling_FisherTest_eventfreq.tsv')
# SSCUS_SNV_100P_TT_df$SNV_freq <- 1
# SSCUS_SNV_10P_TT_df <- Load_TTRes_DF('../TT/Combined Data/SSC/SNV0.1/ SSC_unaffected_sibling_FisherTest_eventfreq_snvfreq0.1.tsv')
# SSCUS_SNV_10P_TT_df$SNV_freq <- 2
# SSCUS_SNV_1P_TT_df <- Load_TTRes_DF('../TT/Combined Data/SSC/SNV0.01/ SSC_unaffected_sibling_FisherTest_eventfreq_snvfreq0.01.tsv')
# SSCUS_SNV_1P_TT_df$SNV_freq <- 3
# 
# SSCUS_TT_df <- rbind(SSCUS_SNV_100P_TT_df,SSCUS_SNV_10P_TT_df,SSCUS_SNV_1P_TT_df)
# SSCUS_text_TT_df <- Get_Interpolated_Text_TT_DF(SSCUS_TT_df, 50)
# 
# Create_TT_Result_Plot(SSCUS_TT_df, SSCUS_text_TT_df, 'SSC Parent-Unaffected Sibling', 'SSC_US_TTP')
