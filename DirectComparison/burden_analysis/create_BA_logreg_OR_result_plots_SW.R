################################################################################################################################################################################
#
# create_BA_logreg_OR_result_plots_SW.R
# purpose: outputs burden analysis logistic regression results as odds ratio and confidence interval plots 
# input: logreg result tables with OR and confidence intervals ending with ...clogit_res.csv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_BA/data/LogRegResultsWithOR
#   -> output of unbiased_burden_analysis_MSSNG+SSC_SW.R 
# output: saves OR vs. SNV freq plots for each data set and comparison to 
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_BA/figures
#   -> Shania final presentation plots
#
################################################################################################################################################################################

library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_BA/data/LogRegResultsWithOR")

Combine_Slide_Data <- function(df.allSNV, df.SNV1p, df.SNV10p){
  df.SNV1p$SNV.freq <- "< 1%"
  df.SNV10p$SNV.freq <- "< 10%"
  df.allSNV$SNV.freq <- "< 100%"
  
  comb.df <- rbind(df.SNV1p, df.SNV10p, df.allSNV)
  return(comb.df)
}

Create_OR_Plot <- function(logreg.df, title, name){
  
  logreg.df <- logreg.df[which(!logreg.df$variant.type %in% c("all_base", "lof")),]
  logreg.df$SNV.freq <- paste("SNV Frequency ", logreg.df$SNV.freq)
  
  ggplot(logreg.df, aes(x = variant.type, y = OR, fill = P < 0.05)) + 
    geom_bar(stat="identity", show.legend = T, color = "black", width = .5) +
    geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), show.legend = F,
                  position = position_dodge(width = .5), size = 1, width = 0.4) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() + ylab("Odds Ratio") + xlab("Variant Type") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 15), 
          axis.text.y = element_text(size = 15), 
          axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          panel.border = element_rect(fill = NA),
          legend.text = element_text(size=13),
          legend.title = element_text(size=13)) +
    coord_cartesian(ylim = c(0, 3)) +
    scale_fill_manual(values = c("lightgrey", "red"), name = "p-value", labels = c(">= 0.05", "< 0.05")) + 
    # scale_fill_viridis(discrete = TRUE, option = "B")+
    scale_x_discrete(labels=c("All", "Missense", "Synonymous")) +
    # scale_y_continuous(breaks = c(1,  5,  9, 15, 22)) +
    # annotate("text", x = c(2.7), y = c(15), angle = 90,
    #          label = paste0("p==", c("2~x~10^-3")), parse = T,
    #          fontface = "italic", cex = 5.5) +
    facet_wrap(.~SNV.freq, nrow=3) +
    theme(strip.text.x = element_text(size = 15, face = "bold.italic"),
          strip.text.y = element_text(size = 15, face = "bold.italic"),
          strip.background = element_rect(color="black", fill="lightgrey", 
                                          size=1.5, linetype="solid")) + 
    ggtitle(sprintf("%s Odds Ratio Plot", title))
  
  ggsave(sprintf("../../figures/%s.OR.plot.png", name), width = 6, height =8)
}
  

## pRec > 0.9 plots
MSSNG.SSC.allSNV.slide.data <- fread("MSSNG.SSC_parent_proband.allSNV.pRec0.9_clogit_res.tsv")
MSSNG.SSC.SNV10p.slide.data <- fread("MSSNG.SSC_parent_proband.SNV10p.pRec0.9_clogit_res.tsv")
MSSNG.SSC.SNV1p.slide.data <- fread("MSSNG.SSC_parent_proband.SNV1p.pRec0.9_clogit_res.tsv")
MSSNG.SSC.slide.data <- Combine_Slide_Data(df.allSNV=MSSNG.SSC.allSNV.slide.data, 
                                            df.SNV1p=MSSNG.SSC.SNV1p.slide.data, 
                                            df.SNV10p=MSSNG.SSC.SNV10p.slide.data)
Create_OR_Plot(MSSNG.SSC.slide.data, 
               title = "MSSNG+SSC pRec > 0.9", name = "MSSNG.SSC.pRec0.9")


SSC_parent_US.allSNV.slide.data <- fread("SSC_parent_US.allSNV.pRec0.9_clogit_res.tsv")
SSC_parent_US.SNV10p.slide.data <-  fread("SSC_parent_US.SNV10p.pRec0.9_clogit_res.tsv")
SSC_parent_US.SNV1p.slide.data <-  fread("SSC_parent_US.SNV1p.pRec0.9_clogit_res.tsv")
SSC_parent_US.slide.data <- Combine_Slide_Data(df.allSNV=SSC_parent_US.allSNV.slide.data,
                                                df.SNV1p = SSC_parent_US.SNV1p.slide.data,
                                                df.SNV10p=SSC_parent_US.SNV10p.slide.data)
Create_OR_Plot(SSC_parent_US.slide.data, 
               title = "SSC parent-US pRec > 0.9", name = "SSC_parent_US.pRec0.9")

SSC_proband_US.allSNV.slide.data <- fread("SSC_proband_US.allSNV.pRec0.9_clogit_res.tsv")
SSC_proband_US.SNV10p.slide.data <-  fread("SSC_proband_US.SNV10p.pRec0.9_clogit_res.tsv")
SSC_proband_US.SNV1p.slide.data <-  fread("SSC_proband_US.SNV1p.pRec0.9_clogit_res.tsv")
SSC_proband_US.slide.data <- Combine_Slide_Data(df.allSNV=SSC_proband_US.allSNV.slide.data,
                                                df.SNV1p = SSC_proband_US.SNV1p.slide.data,
                                                df.SNV10p=SSC_proband_US.SNV10p.slide.data)
Create_OR_Plot(SSC_proband_US.slide.data,
               title = "SSC proband-US pRec > 0.9", name = "SSC_proband_US.pRec0.9")       

## no pRec cut-off plots
MSSNG.SSC.allSNV.nofilt.slide.data <- fread("MSSNG.SSC_parent_proband.allSNV.nofilt_clogit_res.tsv")
MSSNG.SSC.SNV10p.nofilt.slide.data <- fread("MSSNG.SSC_parent_proband.SNV10p.nofilt_clogit_res.tsv")
MSSNG.SSC.SNV1p.nofilt.slide.data <- fread("MSSNG.SSC_parent_proband.SNV1p.nofilt_clogit_res.tsv")
MSSNG.SSC.nofilt.slide.data <- Combine_Slide_Data(df.allSNV=MSSNG.SSC.allSNV.nofilt.slide.data, 
                                           df.SNV1p=MSSNG.SSC.SNV1p.nofilt.slide.data, 
                                           df.SNV10p=MSSNG.SSC.SNV10p.nofilt.slide.data)
Create_OR_Plot(MSSNG.SSC.nofilt.slide.data, 
               title = "MSSNG+SSC parent-proband no pRec filter", name = "MSSNG.SSC")

SSC_parent_US.allSNV.nofilt.slide.data <- fread("SSC_parent_US.allSNV.nofilt_clogit_res.tsv")
SSC_parent_US.SNV10p.nofilt.slide.data <-  fread("SSC_parent_US.SNV10p.nofilt_clogit_res.tsv")
SSC_parent_US.SNV1p.nofilt.slide.data <-  fread("SSC_parent_US.SNV1p.nofilt_clogit_res.tsv")
SSC_parent_US.nofilt.slide.data <- Combine_Slide_Data(df.allSNV=SSC_parent_US.allSNV.nofilt.slide.data,
                                               df.SNV1p = SSC_parent_US.SNV1p.nofilt.slide.data,
                                               df.SNV10p=SSC_parent_US.SNV10p.nofilt.slide.data)
Create_OR_Plot(SSC_parent_US.nofilt.slide.data, 
               title = "SSC parent-US no pRec filter", name = "SSC_parent_US" )

SSC_proband_US.allSNV.nofilt.slide.data <- fread("SSC_proband_US.allSNV.nofilt_clogit_res.tsv")
SSC_proband_US.SNV10p.nofilt.slide.data <-  fread("SSC_proband_US.SNV10p.nofilt_clogit_res.tsv")
SSC_proband_US.SNV1p.nofilt.slide.data <-  fread("SSC_proband_US.SNV1p.nofilt_clogit_res.tsv")
SSC_proband_US.nofilt.slide.data <- Combine_Slide_Data(df.allSNV=SSC_proband_US.allSNV.nofilt.slide.data,
                                                df.SNV1p = SSC_proband_US.SNV1p.nofilt.slide.data,
                                                df.SNV10p=SSC_proband_US.SNV10p.nofilt.slide.data)
Create_OR_Plot(SSC_proband_US.nofilt.slide.data,
               title = "SSC proband-US no pRec filter", name = "SSC_proband_US")     

