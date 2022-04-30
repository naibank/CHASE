library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

#Load("../Unbiased/DT_and_TT.RData")

Load_LogRegRes_DF <- function(path, add_cols = F) {
  var_type <- data.frame('Variant Type'=c('All Variants', 'Synonymous', 'Missense', 'LoF'))
  res <- fread(path, data.table=F)[-5, -1]
  if (add_cols) {
    res$G2 <- res$G2 + res$G3
  }
  res$G1 <- as.character(res$G1)
  res$G2 <- as.character(res$G2)
  res$G3 <- as.character(res$G3)
  return (data.frame(cbind(var_type,res)))
}

Get_Interpolated_Text_DF <- function(df, n=21, skip_lof=F, is_probUS=F) {
  df$P <- -log10(df$P)
  res <- df[, c('SNV_freq', 'P', 'Variant.Type', 'G1', 'G2','G3')]
  res$count_text <- paste('(',res$G1,',',res$G2,',',res$G3,')',sep='')
  if (is_probUS) {
    res$count_text <- paste('(',res$G1,',',res$G2,')',sep='')
  }
  
  df_text_all <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$Variant.Type=='All Variants'], df$P[df$Variant.Type=='All Variants'], n=n)))
  colnames(df_text_all) <- c('SNV_freq', 'P')
  df_text_all$Variant.Type <- 'All Variants'
  df_text_all$G1 <- ""
  df_text_all$G2 <- ""
  df_text_all$G3 <- ""
  df_text_all$count_text <- ""
  
  df_text_syn <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$Variant.Type=='Synonymous'], df$P[df$Variant.Type=='Synonymous'], n=n)))
  colnames(df_text_syn) <- c('SNV_freq', 'P')
  df_text_syn$Variant.Type <- 'Synonymous'
  df_text_syn$G1 <- ""
  df_text_syn$G2 <- ""
  df_text_syn$G3 <- ""
  df_text_syn$count_text <- ""
  
  df_text_mis <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$Variant.Type=='Missense'], df$P[df$Variant.Type=='Missense'], n=n)))
  colnames(df_text_mis) <- c('SNV_freq', 'P')
  df_text_mis$Variant.Type <- 'Missense'
  df_text_mis$G1 <- ""
  df_text_mis$G2 <- ""
  df_text_mis$G3 <- ""
  df_text_mis$count_text <- ""
  
  final <- bind_rows(res, df_text_all, df_text_syn, df_text_mis)
  
  if (!skip_lof) {
    df_text_lof <- as.data.frame(do.call(cbind, approx(df$SNV_freq[df$Variant.Type=='LoF'], df$P[df$Variant.Type=='LoF'], n=n)))
    colnames(df_text_lof) <- c('SNV_freq', 'P')
    df_text_lof$Variant.Type <- 'LoF'
    df_text_lof$G1 <- ""
    df_text_lof$G2 <- ""
    df_text_lof$G3 <- ""
    df_text_lof$count_text <- ""
    final <- bind_rows(final, df_text_lof)
  }
  
  return (final)
}

Create_DT_Result_Plot <- function(df, df_text, title, name, save_path='../DT/FinalPresPlots/') {
  ggplot(data=df, aes(x=SNV_freq, y=-log10(P), color=Variant.Type)) +
    geom_point() +
    # geom_point(data=df[df$P <= 0.05,],
    #            aes(x=SNV_freq,y=coeff,color=Variant.Type),
    #            shape=9,
    #            size=3) +
    geom_line(aes(group=Variant.Type)) +
    scale_x_discrete(limits=c('100%', '10%', '1%')) +
    geom_hline(yintercept=-log10(0.05), linetype='dashed',color='red') +
    geom_text(aes(x=0.55, y=-log10(0.05)+0.05, label='p = 0.05'), color='red', size=2) +
    geom_text_repel(data=df_text,
                    mapping=aes(y=P, x=SNV_freq, label=count_text),
                    force=10,
                    max.overlaps = Inf,
                    size=3,
                    segment.linetype=2) +
    # geom_point(data=df_text, aes(x=SNV_freq, y=P, color=Variant.Type)) +
    labs(y='-log10(P)', x='SNV MAF', color="Variant Type", title=title) + 
    theme(plot.title=element_text(hjust=0.5))
  ggsave(paste(save_path, paste(name,'.png',sep='')), width=6, height=6)
}

#### MSSNG parent-proband Results ####
MSSNG_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/MSSNG_parent_proband_clogit_res.csv')
MSSNG_SNV_100P_regres_df$SNV_freq <- 1
MSSNG_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/MSSNG_parent_proband_SNV10P_clogit_res.csv')
MSSNG_SNV_10P_regres_df$SNV_freq <- 2
MSSNG_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/MSSNG_parent_proband_SNV1P_clogit_res.csv')
MSSNG_SNV_1P_regres_df$SNV_freq <- 3

MSSNG_regres_df <- rbind(MSSNG_SNV_100P_regres_df,MSSNG_SNV_10P_regres_df,MSSNG_SNV_1P_regres_df)
MSSNG_text_regres_df <- Get_Interpolated_Text_DF(MSSNG_regres_df, 50)

# MSSNG_regres_df <- MSSNG_regres_df[MSSNG_regres_df$Variant.Type != 'LoF',]
# MSSNG_text_regres_df <- MSSNG_text_regres_df[MSSNG_text_regres_df$Variant.Type != 'LoF',]

Create_DT_Result_Plot(MSSNG_regres_df, MSSNG_text_regres_df, 'MSSNG Parent-Proband', 'MSSNG_regresP')


#### SSC parent-proband Results ####
SSC_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_parent_proband_clogit_res.csv')
SSC_SNV_100P_regres_df$SNV_freq <- 1
SSC_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_parent_proband_SNV10P_clogit_res.csv')
SSC_SNV_10P_regres_df$SNV_freq <- 2
SSC_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_parent_proband_SNV1P_clogit_res.csv')
SSC_SNV_1P_regres_df$SNV_freq <- 3

SSC_regres_df <- rbind(SSC_SNV_100P_regres_df,SSC_SNV_10P_regres_df,SSC_SNV_1P_regres_df)
SSC_text_regres_df <- Get_Interpolated_Text_DF(SSC_regres_df, 50)

# SSC_regres_df <- SSC_regres_df[SSC_regres_df$Variant.Type != 'LoF',]
# SSC_text_regres_df <- SSC_text_regres_df[SSC_text_regres_df$Variant.Type != 'LoF',]

Create_DT_Result_Plot(SSC_regres_df, SSC_text_regres_df, 'SSC Parent-Proband', 'SSC_regresP')

#### SSC parent-US Results ####
SSCUS_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_parent_US_clogit_res.csv')
SSCUS_SNV_100P_regres_df$SNV_freq <- 1
SSCUS_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_parent_US_SNV10P_clogit_res.csv')
SSCUS_SNV_10P_regres_df$SNV_freq <- 2
SSCUS_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_parent_US_SNV1P_clogit_res.csv')
SSCUS_SNV_1P_regres_df$SNV_freq <- 3

SSCUS_regres_df <- rbind(SSCUS_SNV_100P_regres_df,SSCUS_SNV_10P_regres_df,SSCUS_SNV_1P_regres_df)
SSCUS_text_regres_df <- Get_Interpolated_Text_DF(SSCUS_regres_df, 50)

# SSCUS_regres_df <- SSCUS_regres_df[SSCUS_regres_df$Variant.Type != 'LoF',]
# SSCUS_text_regres_df <- SSCUS_text_regres_df[SSCUS_text_regres_df$Variant.Type != 'LoF',]

Create_DT_Result_Plot(SSCUS_regres_df, SSCUS_text_regres_df, 'SSC Parent-Unaffected Sibling', 'SSC_US_regresP')

#### SSC proband-US Results ####
SSCprobUS_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_proband_US_clogit_res.csv', add_cols = F)
SSCprobUS_SNV_100P_regres_df$SNV_freq <- 1
SSCprobUS_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_proband_US_SNV10P_clogit_res.csv', add_cols = F)
SSCprobUS_SNV_10P_regres_df$SNV_freq <- 2
SSCprobUS_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SSC_proband_US_SNV1P_clogit_res.csv', add_cols = F)
SSCprobUS_SNV_1P_regres_df$SNV_freq <- 3

SSCprobUS_regres_df <- rbind(SSCprobUS_SNV_100P_regres_df,SSCprobUS_SNV_10P_regres_df,SSCprobUS_SNV_1P_regres_df)
SSCprobUS_regres_df <- na.omit(SSCprobUS_regres_df)
SSCprobUS_regres_df <- SSCprobUS_regres_df[!(SSCprobUS_regres_df$G1 == '0' &
                                               SSCprobUS_regres_df$G2 == '0' &
                                               SSCprobUS_regres_df$G3 == '0'),]
# SSCprobUS_regres_df <- SSCprobUS_regres_df[SSCprobUS_regres_df$Variant.Type != 'LoF',]

SSCprobUS_text_regres_df <- Get_Interpolated_Text_DF(SSCprobUS_regres_df, 50, skip_lof=T, is_probUS = T)


Create_DT_Result_Plot(SSCprobUS_regres_df, SSCprobUS_text_regres_df, 'SSC Proband-Unaffected Sibling', 'SSC_probandUS_regresP')


#### SPARK parent-proband Results ####
SPARK_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_parent_proband_clogit_res.csv')
SPARK_SNV_100P_regres_df$SNV_freq <- 1
SPARK_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_parent_proband_SNV10P_clogit_res.csv')
SPARK_SNV_10P_regres_df$SNV_freq <- 2
SPARK_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_parent_proband_SNV1P_clogit_res.csv')
SPARK_SNV_1P_regres_df$SNV_freq <- 3

SPARK_regres_df <- rbind(SPARK_SNV_100P_regres_df,SPARK_SNV_10P_regres_df,SPARK_SNV_1P_regres_df)
SPARK_regres_df <- na.omit(SPARK_regres_df)
SPARK_regres_df <- SPARK_regres_df[!(SPARK_regres_df$G1 == '0' &
                                       SPARK_regres_df$G2 == '0' &
                                       SPARK_regres_df$G3 == '0'),]
SPARK_text_regres_df <- Get_Interpolated_Text_DF(SPARK_regres_df, 50, skip_lof = T)

# SPARK_regres_df <- SPARK_regres_df[SPARK_regres_df$Variant.Type != 'LoF',]
# SPARK_text_regres_df <- SPARK_text_regres_df[SPARK_text_regres_df$Variant.Type != 'LoF',]

Create_DT_Result_Plot(SPARK_regres_df, SPARK_text_regres_df, 'SPARK Parent-Proband', 'SPARK_regresP')

#### SPARK parent-US Results ####
SPARKUS_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_parent_US_clogit_res.csv')
SPARKUS_SNV_100P_regres_df$SNV_freq <- 1
SPARKUS_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_parent_US_SNV10P_clogit_res.csv')
SPARKUS_SNV_10P_regres_df$SNV_freq <- 2
SPARKUS_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_parent_US_SNV1P_clogit_res.csv')
SPARKUS_SNV_1P_regres_df$SNV_freq <- 3

SPARKUS_regres_df <- rbind(SPARKUS_SNV_100P_regres_df,SPARKUS_SNV_10P_regres_df,SPARKUS_SNV_1P_regres_df)
SPARKUS_regres_df <- na.omit(SPARKUS_regres_df)
SPARKUS_regres_df <- SPARKUS_regres_df[!(SPARKUS_regres_df$G1 == '0' &
                                           SPARKUS_regres_df$G2 == '0' &
                                           SPARKUS_regres_df$G3 == '0'),]
SPARKUS_text_regres_df <- Get_Interpolated_Text_DF(SPARKUS_regres_df, 50, skip_lof = T)

# SPARKUS_regres_df <- SPARKUS_regres_df[SPARKUS_regres_df$Variant.Type != 'LoF',]
# SPARKUS_text_regres_df <- SPARKUS_text_regres_df[SPARKUS_text_regres_df$Variant.Type != 'LoF',]

Create_DT_Result_Plot(SPARKUS_regres_df, SPARKUS_text_regres_df, 'SPARK Parent-Unaffected Sibling', 'SPARK_US_regresP')

#### SPARK proband-US Results ####
SPARKprobUS_SNV_100P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_proband_US_clogit_res.csv', add_cols = F)
SPARKprobUS_SNV_100P_regres_df$SNV_freq <- 1
SPARKprobUS_SNV_10P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_proband_US_SNV10P_clogit_res.csv', add_cols = F)
SPARKprobUS_SNV_10P_regres_df$SNV_freq <- 2
SPARKprobUS_SNV_1P_regres_df <- Load_LogRegRes_DF('../DT/LogRegResults/SPARK_proband_US_SNV1P_clogit_res.csv', add_cols = F)
SPARKprobUS_SNV_1P_regres_df$SNV_freq <- 3

SPARKprobUS_regres_df <- rbind(SPARKprobUS_SNV_100P_regres_df,SPARKprobUS_SNV_10P_regres_df,SPARKprobUS_SNV_1P_regres_df)
SPARKprobUS_regres_df <- na.omit(SPARKprobUS_regres_df)
SPARKprobUS_regres_df <- SPARKprobUS_regres_df[!(SPARKprobUS_regres_df$G1 == '0' &
                                               SPARKprobUS_regres_df$G2 == '0' &
                                               SPARKprobUS_regres_df$G3 == '0'),]
# SPARKprobUS_regres_df <- SPARKprobUS_regres_df[SPARKprobUS_regres_df$Variant.Type != 'LoF',]

SPARKprobUS_text_regres_df <- Get_Interpolated_Text_DF(SPARKprobUS_regres_df, 50, skip_lof=T, is_probUS = T)


Create_DT_Result_Plot(SPARKprobUS_regres_df, SPARKprobUS_text_regres_df, 'SPARK Proband-Unaffected Sibling', 'SPARK_probandUS_regresP')
