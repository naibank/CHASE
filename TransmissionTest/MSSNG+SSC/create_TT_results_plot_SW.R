################################################################################################################################################################################
# 
# create_TT_results_plot_SW.R
# purpose: plots MSSNG+SSC transmission Fisher's exact test results as var. type vs. transmission rate
# input: Fisher's test results, e.g. MSSNG.SSC_parent_proband_FisherTest_snvfreq1.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data
#             -> output of MSSNG+SSC_get_FisherTest_snvfreq_SW.R
# output: variant type vs. transmission rates for SNV frequency cut-offs of 100%, 10%, 1%
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/figures
#             -> Shania final pres. plots
# 
# notes: plots save as "OR plots" due to a previous version, but are actually "transmission rate plots"
#
##############################################################################################################################################################################

library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data")

Load_TTRes_DF <- function(path, add_cols = F, pRec.filter = T) {
  res <- fread(path, data.table=F)
  if (pRec.filter){
    res <- res[res$pRec == 0.9,] 
  }
  if (!pRec.filter){
    res <- res[res$pRec == 0.0,] 
  }
  res$target_UT <- res$target.parent - res$target.child
  res$target.child <- as.character(res$target.child)
  res$target_UT <- as.character(res$target_UT)
  
  return (data.frame(res))
}

Combine_Slide_Data <- function(df.allSNV, df.SNV1p, df.SNV10p){
  df.SNV1p$SNV.freq <- "< 1%"
  df.SNV10p$SNV.freq <- "< 10%"
  df.allSNV$SNV.freq <- "< 100%"
  
  comb.df <- rbind(df.SNV1p, df.SNV10p, df.allSNV)
  return(comb.df)
}

Create_OR_Plot <- function(fisher.df, title, name){
  
  fisher.df$SNV.freq <- paste("SNV Frequency ", fisher.df$SNV.freq)
  ggplot(fisher.df, aes(x = factor(variant_type, c("LoF", "damaging_missense", "nonsynonymous", "synonymous")),
                        y = rate, fill = (P < 0.05 & rate > 0.5))) +
    geom_bar(stat="identity", show.legend = T, color = "black", width = .5) +
    geom_errorbar(aes(ymin = rate_lower, ymax = rate_upper), show.legend = F,
                  position = position_dodge(width = .5), size = 1, width = 0.4) +
    geom_hline(yintercept = 0.5, lty = 2) +
    theme_classic() + ylab("Transmission rate") + xlab("Variant Type") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 15), 
          axis.text.y = element_text(size = 15), 
          axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          panel.border = element_rect(fill = NA),
          legend.text = element_text(size=13),
          legend.title = element_text(size=13)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = c("lightgrey", "red"), name = "p-value", labels = c(">= 0.05", "< 0.05")) + 
    scale_x_discrete(labels=c("LoF", "Damaging Missense", "Nonsynonymous", "Synonymous")) +
    facet_wrap(.~SNV.freq, nrow=3) +
    theme(strip.text.x = element_text(size = 15, face = "bold.italic"),
          strip.text.y = element_text(size = 15, face = "bold.italic"),
          strip.background = element_rect(color="black", fill="lightgrey", 
                                          size=1.5, linetype="solid")) + 
    ggtitle(sprintf("%s Odds Ratio Plot", title))
  
  ggsave(sprintf("../figures/%s.OR.plot.png", name), width = 6, height =8)
}


##############################################################################################################################################################################
##### MSSNG parent-proband results ####
## pRec > 0.9
MSSNG.SSC_SNV_100P_TT_df <- Load_TTRes_DF('./MSSNG.SSC_parent_proband_FisherTest_snvfreq1.tsv')
MSSNG.SSC_SNV_10P_TT_df <- Load_TTRes_DF('./MSSNG.SSC_parent_proband_FisherTest_snvfreq0.1.tsv')
MSSNG.SSC_SNV_1P_TT_df <- Load_TTRes_DF('./MSSNG.SSC_parent_proband_FisherTest_snvfreq0.01.tsv')

MSSNG.SSC_TT_pRec0.9_df <- Combine_Slide_Data(df.allSNV=MSSNG.SSC_SNV_100P_TT_df, 
                                              df.SNV1p=MSSNG.SSC_SNV_1P_TT_df, 
                                              df.SNV10p=MSSNG.SSC_SNV_10P_TT_df)
Create_OR_Plot(MSSNG.SSC_TT_pRec0.9_df,
               title = "MSSNG+SSC Parent-Proband FisherTest pRec > 0.9", 
               name = "MSSNG.SSC_parent_proband_FisherTest.pRec") 

## no pRec cut-off
MSSNG.SSC_SNV_100P_TT_nopRec_df <- Load_TTRes_DF('./MSSNG.SSC_parent_proband_FisherTest_snvfreq1.tsv', pRec.filter = F)
MSSNG.SSC_SNV_10P_TT_nopRec_df <- Load_TTRes_DF('./MSSNG.SSC_parent_proband_FisherTest_snvfreq0.1.tsv', pRec.filter = F)
MSSNG.SSC_SNV_1P_TT_nopRec_df <- Load_TTRes_DF('./MSSNG.SSC_parent_proband_FisherTest_snvfreq0.01.tsv', pRec.filter = F)

MSSNG.SSC_TT_nopRec_df <- Combine_Slide_Data(df.allSNV=MSSNG.SSC_SNV_100P_TT_nopRec_df, 
                                             df.SNV1p=MSSNG.SSC_SNV_1P_TT_nopRec_df, 
                                             df.SNV10p=MSSNG.SSC_SNV_10P_TT_nopRec_df)
Create_OR_Plot(MSSNG.SSC_TT_nopRec_df,
               title = "MSSNG+SSC Parent-Proband FisherTest all pRec", 
               name = "MSSNG.SSC_parent_proband_FisherTest.nopRec")  


##############################################################################################################################################################################
#### SSC parent-US Results ####
## pRec > 0.9
SSCUS_SNV_100P_TT_df <- Load_TTRes_DF('./SSC_parent_US_FisherTest_snvfreq1.tsv')
SSCUS_SNV_10P_TT_df <- Load_TTRes_DF('./SSC_parent_US_FisherTest_snvfreq0.1.tsv')
SSCUS_SNV_1P_TT_df <- Load_TTRes_DF('./SSC_parent_US_FisherTest_snvfreq0.01.tsv')

SSCUS_TT_pRec0.9_df <- Combine_Slide_Data(df.allSNV=SSCUS_SNV_100P_TT_df,
                                          df.SNV1p=SSCUS_SNV_1P_TT_df,
                                          df.SNV10p=SSCUS_SNV_10P_TT_df)
Create_OR_Plot(SSCUS_TT_pRec0.9_df,
               title = "SSC Parent-US FisherTest pRec > 0.9", 
               name = "SSC_praent_unaffectedsib_FisherTest.pRec") 

## no pRec cut-off
SSCUS_SNV_100P_TT_nopRec_df <- Load_TTRes_DF('./SSC_parent_US_FisherTest_snvfreq1.tsv', pRec.filter = F)
SSCUS_SNV_10P_TT_nopRec_df <- Load_TTRes_DF('./SSC_parent_US_FisherTest_snvfreq0.1.tsv', pRec.filter = F)
SSCUS_SNV_1P_TT_nopRec_df <- Load_TTRes_DF('./SSC_parent_US_FisherTest_snvfreq0.01.tsv', pRec.filter = F)

SSCUS_TT_nopRec_df <- Combine_Slide_Data(df.allSNV=SSCUS_SNV_100P_TT_nopRec_df,
                                          df.SNV1p=SSCUS_SNV_1P_TT_nopRec_df,
                                          df.SNV10p=SSCUS_SNV_10P_TT_nopRec_df)
Create_OR_Plot(SSCUS_TT_nopRec_df,
               title = "SSC Parent-US FisherTest all pRec", 
               name = "SSC_praent_unaffectedsib_FisherTest.nopRec") 
