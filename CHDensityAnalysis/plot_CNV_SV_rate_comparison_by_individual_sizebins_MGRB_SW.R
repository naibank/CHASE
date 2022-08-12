################################################################################################################################################################################
# 
# plot_CNV_SV_rate_comparison_by_individual_sizebins_MGRB.R
# purpose: plots MGRB sliding window cut-offs (median, 80, 90 perc.) spline
#        /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/CH_hits/
# output: sliding.window.sizebins.MGRB_median.mean.80.90.individuals_2.5kb.png
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/figures/
# 
# notes: 
#
##############################################################################################################################################################################


library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggformula)
library(broom)
library(cowplot)
library(reshape2)

setwd("/Users/shaniawu/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette.t <- c("#99999933", "#E69F0033", "#56B4E933", "#009E7333", "#F0E44233", "#0072B233", "#D55E0033", "#CC79A733")
cbVermeer <- c('#6495ED','#93CCEA','#6495ED50','#93CCEA50','#6495ED20','#93CCEA20')

## Read CH hits table
MGRB.CH.hits <- fread("./data/CH_hits/MGRB_CH_hits_nofilt.tsv", data.table=F)


########################################################################################################################################################################################
#### Get CUT-OFFs for each size bin ####

Get_CutOffs_Table <- function(CH.hits, bin.end=F, bin.size){
  ## remove (0,0) individuals from CH_hits (outliers)
  CH_hits <- CH_hits[which(!(CH_hits$exSize == 0 & CH_hits$CH_hit == 0)),]
  
  ## make size.bins list
  size.bins <- list()
  bin.start <- 0
  if (bin.end == F){
    bin.end <- max(CH_hits$exSize)
  }
  while (bin.start < bin.end){
    bin <- list(bin.start, bin.start + bin.size) # list(start, end)
    size.bins <- c(size.bins, list(bin))
    bin.start <- bin.start + bin.size
  }
  
  ## find cut-offs for each size bin
  cutoffs.out <- data.frame()
  
  for (i in 1:length(size.bins)){
    bin.start.end <- size.bins[[i]]
    start <- bin.start.end[[1]]
    end <- bin.start.end[[2]]
    CH_hits_in_bin <- CH_hits[which(CH_hits$exSize >= start & CH_hits$exSize < end),] 
    
    ## add size bin start & end columns
    cutoffs.df <- data.frame(size.bin.start=start)
    cutoffs.df$size.bin.end <- end
    
    ## add cutoff columns
    cutoffs.df$median <- median(CH_hits_in_bin$CH_hit)
    cutoffs.df$mean <- mean(CH_hits_in_bin$CH_hit)
    cutoffs.df$percentile_0.8 <- quantile(CH_hits_in_bin$CH_hit, probs = 0.8)
    cutoffs.df$percentile_0.9 <- quantile(CH_hits_in_bin$CH_hit, probs = 0.9)
    cutoffs.df$bin.count <- nrow(CH_hits_in_bin)
 
    cutoffs.out <- rbind(cutoffs.out, cutoffs.df)
  }
  return(cutoffs.out)
}

MGRB.cutoffs.df <- Get_CutOffs_Table(MGRB.CH.hits, bin.end = 30000, bin.size = 2500)


#### Plot CUT-OFFs for each size bin - sliding window approach ####

Get_Sliding_Window_Cutoff_Plots <- function(cutoffs.df, CH.hits.df, file.name, plot.individuals=F){
  ### Plots median, mean, 80th, and 90th percentiles for each 2.5kb size bin 
  ###   using sliding window approach (takes average of 3 size bins)
  
  ## Make sliding window df with window averages
  sliding.window.df <- data.frame()
  
  for (i in 2:(nrow(cutoffs.df)-3)){
    if (cutoffs.df[i,]$size.bin.start == 20000){ # get average of all size bins >20kb & <30kb
      window.start.i <- i 
      window.end.i <- i+3 
      cutoffs.df[i,]$size.bin.end <- 30000
    } else {
      window.start.i <- i-1
      window.end.i <- i+1
    }
    window.median <- mean(cutoffs.df$median[window.start.i:window.end.i], na.rm=TRUE)
    window.mean <- mean(cutoffs.df$mean[window.start.i:window.end.i], na.rm=TRUE)
    window.80th <- mean(cutoffs.df$"percentile_0.8"[window.start.i:window.end.i],na.rm=TRUE)
    window.90th <- mean(cutoffs.df$"percentile_0.9"[window.start.i:window.end.i],na.rm=TRUE)
    sliding.window.df.row <- cutoffs.df[i,]
    sliding.window.df.row$median <- window.median
    sliding.window.df.row$mean <- window.mean
    sliding.window.df.row$`percentile_0.8` <- window.80th
    sliding.window.df.row$`percentile_0.9` <- window.90th
    sliding.window.df <- rbind(sliding.window.df, sliding.window.df.row)
  }
  
  CH.hits.df <- CH.hits.df[which(CH.hits.df$exSize <= 30000),]
      
  ## Plot sliding window averages
  sliding.window.df.melt <- melt(sliding.window.df, 
                                 id.vars = c("size.bin.start","size.bin.end", "bin.count"))
  ggplot(sliding.window.df.melt, #[which(sliding.window.df.melt$variable == "median"),]
         aes(x = size.bin.start, y = value)) +
    geom_point(aes(colour = variable), size = 2) + 
    geom_point(data=CH.hits.df, aes(x = exSize, y = CH_hit, alpha = 0.1)) +
    scale_color_manual(values = c("percentile_0.9"= "#b79f01",
                                  "percentile_0.8"= "#f8766d",
                                  "mean"= "#00ba38",
                                  "median"="#C77CFF")) +
    labs(y='Control No. SNV Varaible', x='Size Bin Start (bp)',
         title="MGRB Sliding Window (bin size = 2.5 kb)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    scale_x_continuous(breaks= seq(0,30000,2500)) +
    geom_spline(data =sliding.window.df.melt, aes(colour = variable))
    # 
    # geom_smooth(data = sliding.window.df.melt, method='lm', alpha = 0.1, aes(colour = variable)) +
    # geom_smooth(data = sliding.window.df.melt, method='lm', alpha = 0.1, linetype = "dashed",
    #             aes(colour = variable), se = FALSE, fullrange = TRUE) 
    # 
  plot.path <- sprintf("./figures/sliding.window.sizebins.%s_2.5kb.png", file.name)
  ggsave(plot.path, width = 9, height = 7)
}


# Get_Sliding_Window_Cutoff_Plots(MSSNG.SSC.median.80.90, "median.80.90")
# Get_Sliding_Window_Cutoff_Plots(MGRB.cutoffs.df, MGRB.CH.hits, "MGRB_median.mean.80.90.spline")
Get_Sliding_Window_Cutoff_Plots(MGRB.cutoffs.df, MGRB.CH.hits, "MGRB_median.mean.80.90.spline.individuals")
