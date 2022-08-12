################################################################################################################################################################################
# 
# create_association_results_plot_SW.R
# purpose: plots MSSNG+SSC spline association analysis Fisher's test results 
# input: Fisher's test results, e.g. MSSNG.SSC.fisher.spline.allcutoffs.pRec.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/
# output: variant type vs. OR plots for different cut-offs: median, 80th & 90th percentiles
#         ../figures/
# 
# notes: SNV freq not tested due to limited no. samples
#
##############################################################################################################################################################################

library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/")


Create_OR_Plot <- function(fisher.df, title, name){

  ggplot(fisher.df, aes(x = factor(cut.off, c("median", "percentile_0.8", "percentile_0.9")),
                        y = OR, fill = p.value < 0.05)) +
    geom_bar(stat="identity", show.legend = T, color = "black", width = .5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), show.legend = F,
                  position = position_dodge(width = .5), size = 1, width = 0.4) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() + ylab("Odds Ratio") + xlab("Cut-off for Fisher Test") +
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
    scale_x_discrete(labels=c("median", "80th percentile", "90th percentile")) +
    # scale_y_continuous(breaks = c(1,  5,  9, 15, 22)) +
    # annotate("text", x = c(2.7), y = c(15), angle = 90,
    #          label = paste0("p==", c("2~x~10^-3")), parse = T,
    #          fontface = "italic", cex = 5.5) +
    facet_wrap(.~deletion.size.range, nrow=2) +
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
assoc.fisher.pRec <- fread("MSSNG.SSC.fisher.spline.allcutoffs.pRec.tsv", data.table=F)

Create_OR_Plot(assoc.fisher.pRec,
               title = "MSSNG+SSC Association Analysis Spline - pRec > 0.9", 
               name = "MSSNG.SSC.assoc.spline.FisherTest.pRec") 

## no pRec cut-off
assoc.fisher.nopRec <- fread("MSSNG.SSC.fisher.spline.allcutoffs.tsv", data.table=F)
assoc.fisher.nopRec$cut.off[which(assoc.fisher.nopRec$cut.off == "80th percentile")] <- "percentile_0.8"
assoc.fisher.nopRec$cut.off[which(assoc.fisher.nopRec$cut.off == "90th percentile")] <- "percentile_0.9"

Create_OR_Plot(assoc.fisher.nopRec,
               title = "MSSNG+SSC Association Analysis Spline - all pRec", 
               name = "MSSNG.SSC.assoc.spline.FisherTest.nopRec")  


