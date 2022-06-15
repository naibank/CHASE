library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)
library(cowplot)

setwd("/Users/shaniawu/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size")

########################################################################################################################################################################################
#### Plot CUT-OFFs for each size bin ####

## median, bin size = 1kb
# MSSNG.SSC <- fread("./data/median.fisher/MSSNG.SSC.median.fisher.tsv", data.table=F)
# MSSNG.SSC.less.30k <- MSSNG.SSC[which(MSSNG.SSC$size.bin.end <= 30000),]
# 
# MSSNG <- fread("./data/median.fisher/MSSNG.median.fisher.tsv", data.table=F)
# MSSNG.less.30k <- MSSNG[which(MSSNG$size.bin.end <= 30000),]
# 
# SSC <- fread("./data/median.fisher/SSC.median.fisher.tsv", data.table=F)
# SSC.less.30k <- SSC[which(SSC$size.bin.end <= 30000),]


# all.files <- list.files(path = folder, pattern = "^filename\\d+_\\d{4}\\.csv$")
# 
# l <- lapply(all.files, fread, sep=",")
# dt <- rbindlist(l )
# setkey(dt, ID)
# unique(dt, by = 'ID')

## median, bin size = 2.5kb
MSSNG.SSC.median <- fread("./data/sizebins.fisher/MSSNG.SSC.median.2.5kb.fisher.tsv", data.table=F)
MSSNG.SSC.median$cohort <- "MSSNG+SSC"
MSSNG.median <- fread("./data/sizebins.fisher/MSSNG.median.2.5kb.fisher.tsv", data.table=F)
MSSNG.median$cohort <- "MSSNG"
SSC.median <- fread("./data/sizebins.fisher/SSC.median.2.5kb.fisher.tsv", data.table=F)
SSC.median$cohort <- "SSC"

comb.median <- rbind(MSSNG.SSC.median, MSSNG.median, SSC.median)
comb.median <- comb.median[which(comb.median$size.bin.end <= 30000),]

ggplot(comb.median, aes(x = size.bin.start, y = cut.off)) +
  geom_point(aes(color=cohort)) +
  geom_line(aes(color=cohort)) +
  scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) +
  labs(y='Control No. SNV Meidan', x='Size Bin Start (bp)', color="cohort", 
       title="MSSNG+SSC Bin Median (bin size = 2.5 kb)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

plot.path <- "./figures/sizebins.median_2.5kb.png"
ggsave(plot.path, width = 6, height = 5)


## mean, bin size = 2.5kb
MSSNG.SSC.mean <- fread("./data/sizebins.fisher/MSSNG.SSC.mean.2.5kb.fisher.tsv", data.table=F)
MSSNG.SSC.mean$cohort <- "MSSNG+SSC"
MSSNG.mean <- fread("./data/sizebins.fisher/MSSNG.mean.2.5kb.fisher.tsv", data.table=F)
MSSNG.mean$cohort <- "MSSNG"
SSC.mean <- fread("./data/sizebins.fisher/SSC.mean.2.5kb.fisher.tsv", data.table=F)
SSC.mean$cohort <- "SSC"

comb.mean <- rbind(MSSNG.SSC.mean, MSSNG.mean, SSC.mean)
comb.mean <- comb.mean[which(comb.mean$size.bin.end <= 30000),]

ggplot(comb.mean, aes(x = size.bin.start, y = cut.off)) +
  geom_point(aes(color=cohort)) +
  geom_line(aes(color=cohort)) +
  scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) +
  labs(y='Control No. SNV Mean', x='Size Bin Start (bp)', color="cohort", 
       title="MSSNG+SSC Bin Mean (bin size = 2.5 kb)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

plot.path <- "./figures/sizebins.mean_2.5kb.png"
ggsave(plot.path, width = 6, height = 5)


## 70th Percentile, bin size = 2.5kb
MSSNG.SSC.0.7 <- fread("./data/sizebins.fisher/MSSNG.SSC.0.7_2.5kb.fisher.tsv", data.table=F)
MSSNG.SSC.0.7$cohort <- "MSSNG+SSC"
MSSNG.0.7 <- fread("./data/sizebins.fisher/MSSNG.0.7_2.5kb.fisher.tsv", data.table=F)
MSSNG.0.7$cohort <- "MSSNG"
SSC.0.7 <- fread("./data/sizebins.fisher/SSC.0.7_2.5kb.fisher.tsv", data.table=F)
SSC.0.7$cohort <- "SSC"

comb.0.7 <- rbind(MSSNG.SSC.0.7, MSSNG.0.7, SSC.0.7)
comb.0.7 <- comb.0.7[which(comb.0.7$size.bin.end <= 30000),]

ggplot(comb.0.7, aes(x = size.bin.start, y = cut.off)) +
  geom_point(aes(color=cohort)) +
  geom_line(aes(color=cohort)) +
  scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) +
  labs(y='Control No. SNV 70th Percentile', x='Size Bin Start (bp)', color="cohort", 
       title="MSSNG+SSC Bin 70th Percentile (bin size = 2.5 kb)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

plot.path <- "./figures/sizebins.0.7_2.5kb.png"
ggsave(plot.path, width = 6, height = 5)


## 80th Percentile, bin size = 2.5kb
MSSNG.SSC.0.8 <- fread("./data/sizebins.fisher/MSSNG.SSC.0.8_2.5kb.fisher.tsv", data.table=F)
MSSNG.SSC.0.8$cohort <- "MSSNG+SSC"
MSSNG.0.8 <- fread("./data/sizebins.fisher/MSSNG.0.8_2.5kb.fisher.tsv", data.table=F)
MSSNG.0.8$cohort <- "MSSNG"
SSC.0.8 <- fread("./data/sizebins.fisher/SSC.0.8_2.5kb.fisher.tsv", data.table=F)
SSC.0.8$cohort <- "SSC"

comb.0.8 <- rbind(MSSNG.SSC.0.8, MSSNG.0.8, SSC.0.8)
comb.0.8 <- comb.0.8[which(comb.0.8$size.bin.end <= 30000),]

ggplot(comb.0.8, aes(x = size.bin.start, y = cut.off)) +
  geom_point(aes(color=cohort)) +
  geom_line(aes(color=cohort)) +
  scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) +
  labs(y='Control No. SNV 80th Percentile', x='Size Bin Start (bp)', color="cohort", 
       title="MSSNG+SSC Bin 80th Percentile (bin size = 2.5 kb)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

plot.path <- "./figures/sizebins.0.8_2.5kb.png"
ggsave(plot.path, width = 6, height = 5)

## 90th Percentile, bin size = 2.5kb
MSSNG.SSC.0.9 <- fread("./data/sizebins.fisher/MSSNG.SSC.0.9_2.5kb.fisher.tsv", data.table=F)
MSSNG.SSC.0.9$cohort <- "MSSNG+SSC"
MSSNG.0.9 <- fread("./data/sizebins.fisher/MSSNG.0.9_2.5kb.fisher.tsv", data.table=F)
MSSNG.0.9$cohort <- "MSSNG"
SSC.0.9 <- fread("./data/sizebins.fisher/SSC.0.9_2.5kb.fisher.tsv", data.table=F)
SSC.0.9$cohort <- "SSC"

comb.0.9 <- rbind(MSSNG.SSC.0.9, MSSNG.0.9, SSC.0.9)
comb.0.9 <- comb.0.9[which(comb.0.9$size.bin.end <= 30000),]

ggplot(comb.0.9, aes(x = size.bin.start, y = cut.off)) +
  geom_point(aes(color=cohort)) +
  geom_line(aes(color=cohort)) +
  scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) +
  labs(y='Control No. SNV 90th Percentile', x='Size Bin Start (bp)', color="cohort", 
       title="MSSNG+SSC Bin 90th Percentile (bin size = 2.5 kb)") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

plot.path <- "./figures/sizebins.0.9_2.5kb.png"
ggsave(plot.path, width = 6, height = 5)



########################################################################################################################################################################################
#### Plot HISTOGRAM & VIOLIN plots ####
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette.t <- c("#99999933", "#E69F0033", "#56B4E933", "#009E7333", "#F0E44233", "#0072B233", "#D55E0033", "#CC79A733")
cbVermeer <- c('#6495ED','#93CCEA','#6495ED50','#93CCEA50','#6495ED20','#93CCEA20')


MSSNG.SSC.CH.hits <- fread("./data/CH_hits/MSSNG_SSC_CH_hits_nofilt.tsv", data.table=F)
MSSNG.SSC.CH.hits$cohort <- "MSSNG+SSC"

MSSNG.CH.hits <- fread("./data/CH_hits/MSSNG_CH_hits_nofilt.tsv", data.table=F)
MSSNG.CH.hits$cohort <- "MSSNG"

SSC.CH.hits <- fread("./data/CH_hits/SSC_CH_hits_nofilt.tsv", data.table=F)
SSC.CH.hits$cohort <- "SSC"

comb.CH.hits <- rbind(MSSNG.SSC.CH.hits, MSSNG.CH.hits, SSC.CH.hits)


Get_Histogram_Violin_Plots <- function(CH_hits_df, size.bin){
  CH_hits_df <- CH_hits_df[which(!(CH_hits_df$exSize == 0 & CH_hits_df$CH_hit == 0)),]
  
  ### histogram
  if (size.bin == "less.5kb"){
    his.plot <- ggplot(CH_hits_df, aes(x = CH_hit)) + 
      geom_histogram(binwidth = 1, aes(fill=cohort), position="dodge", color="black") + 
      coord_cartesian(xlim=c(0,10)) +
      labs(y='No. Individuals', x="No. CH hits")+
      theme(axis.title = element_text(size = 12)) +
      scale_x_continuous(breaks= c(0:20)) 
      
  } else {
    his.plot <- ggplot(CH_hits_df, aes(x = CH_hit)) + 
      geom_histogram(binwidth = 1, aes(fill=cohort), position="dodge", color="black") + 
      coord_cartesian(xlim=c(0,12)) +
      labs(y='No. Individuals', x="No. CH hits")+
      theme(axis.title = element_text(size = 12)) +
      scale_x_continuous(breaks= c(0:12)) 
  }
  ### boxplot + violin plot
  violin <- ggplot(CH_hits_df, aes(y = CH_hit, x = cohort)) + 
    geom_violin(aes(fill=cohort), position="dodge", color="black") +
    geom_boxplot(fill="white", position="dodge", color="black", width = .2, alpha=.5) +
    labs(y='No. CH hits') +
    theme(axis.title = element_text(size = 12))  
  
  ## combine plots
  plot_grid(his.plot, violin, nrow= 1)
  plot.path <- sprintf("./figures/histograms/MSSNG.SSC.hist.violin.%s.png", size.bin)
  ggsave(plot.path, width = 10, height = 5)
  
}

Get_Mean_Var_DI <- function(CH_hits_df, size.bin){
  CH_hits_df <- CH_hits_df[which(!(CH_hits_df$exSize == 0 & CH_hits_df$CH_hit == 0)),]
  
  table <- data.frame()
  
  for (x in c("MSSNG+SSC", "MSSNG", "SSC")){
    CH_hits_df.cohort <- CH_hits_df[which((CH_hits_df$cohort == x)),]
    mean <- signif(mean(CH_hits_df.cohort$CH_hit),4)
    variance <- signif(var(CH_hits_df.cohort$CH_hit),4)
    DI <- signif(variance/mean, 4)
    row <- data.table("Cohort" = x, "Mean" = mean, "Variance" = variance, "Dispersion Index"= DI)
    table <- rbind(table, row)
  }
  
  file.path <- sprintf("./data/stats.%s.tsv", size.bin)
  write.table(table, file.path, sep="\t", row.names=F, quote=F, col.names=T)
}


comb.CH.hits.less.5kb <- comb.CH.hits[which(comb.CH.hits$exSize <= 5000), ]
Get_Histogram_Violin_Plots(comb.CH.hits.less.5kb, "less.5kb")
Get_Mean_Var_DI(comb.CH.hits.less.5kb, "less.5kb")

comb.CH.hits.more.5kb.less.10kb <- comb.CH.hits[which(comb.CH.hits$exSize > 5000 & comb.CH.hits$exSize <= 10000), ]
Get_Histogram_Violin_Plots(comb.CH.hits.more.5kb.less.10kb, "more.5kb.less.10kb")
Get_Mean_Var_DI(comb.CH.hits.more.5kb.less.10kb, "more.5kb.less.10kb")

comb.CH.hits.more.10kb.less.15kb <- comb.CH.hits[which(comb.CH.hits$exSize > 10000 & comb.CH.hits$exSize <= 15000), ]
Get_Histogram_Violin_Plots(comb.CH.hits.more.10kb.less.15kb, "more.10kb.less.15kb")
Get_Mean_Var_DI(comb.CH.hits.more.10kb.less.15kb, "more.10kb.less.15kb")

comb.CH.hits.more.15kb.less.20kb <- comb.CH.hits[which(comb.CH.hits$exSize > 15000 & comb.CH.hits$exSize <= 20000), ]
Get_Histogram_Violin_Plots(comb.CH.hits.more.15kb.less.20kb, "more.15kb.less.20kb")
Get_Mean_Var_DI(comb.CH.hits.more.15kb.less.20kb, "more.15kb.less.20kb")

comb.CH.hits.more.20kb.less.25kb <- comb.CH.hits[which(comb.CH.hits$exSize > 20000 & comb.CH.hits$exSize <= 25000), ]
Get_Histogram_Violin_Plots(comb.CH.hits.more.20kb.less.25kb, "more.20kb.less.25kb")
Get_Mean_Var_DI(comb.CH.hits.more.20kb.less.25kb, "more.20kb.less.25kb")


### Plot Medians (cut-off) ####
# MSSNG.SSC.median.plot <- ggplot(MSSNG.SSC.less.30k, aes(x = size.bin.start, y = cut.off)) +
#   geom_point() +
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#                      axis.title.x = element_text(size = 10),
#                      plot.title = element_text(face = "bold", size = 11)) +
#   scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) + 
#   ggtitle("MSSNG+SSC Bin Median") +
#   xlab("Size Bin Start (bp)") + ylab("Control Median")
# 
# MSSNG.median.plot <- ggplot(MSSNG.less.30k, aes(x = size.bin.start, y = cut.off)) +
#   geom_point() +
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#                      axis.title.x = element_text(size = 10),
#                      plot.title = element_text(face = "bold", size = 11)) +
#   scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) + 
#   ggtitle("MSSNG Bin Median") +
#   xlab("Size Bin Start (bp)") + ylab("Control Median")
# 
# SSC.median.plot <- ggplot(SSC.less.30k, aes(x = size.bin.start, y = cut.off)) +
#   geom_point() +
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#                      axis.title.x = element_text(size = 10),
#                      plot.title = element_text(face = "bold", size = 11)) +
#   scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) + 
#   ggtitle("SSC Bin Median") +
#   xlab("Size Bin Start (bp)") + ylab("Control Median")
# 
# ## Combine plots
# plot_grid(MSSNG.SSC.median.plot, MSSNG.median.plot, SSC.median.plot, nrow= 3)
# plot.path <- "./figures/MSSNG.SSC.sizebins.medianplot.png"
# ggsave(plot.path, width = 6, height = 7)
# 
# 
# 
# ############################################################################################
# ### Plot Odds Ratio ####
# 
# ## MSSNG + SSC
# MSSNG.SSC.plot <- ggplot(MSSNG.SSC.less.30k, aes(x = size.bin.start, y = estimate)) +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.3))+ 
#   coord_cartesian(ylim=c(0,5)) +
#   scale_alpha(guide = 'none') +
#   geom_point() +
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#                      axis.title.x = element_text(size = 10),
#                      plot.title = element_text(face = "bold", size = 11)) +
#   scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) + 
#   ggtitle("MSSNG+SSC") +
#   xlab("Size Bin Start (bp)") + ylab("OR")
# 
# MSSNG.plot <- ggplot(MSSNG.less.30k, aes(x = size.bin.start, y = estimate)) +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.3))+ 
#   coord_cartesian(ylim=c(0,5)) +
#   scale_alpha(guide = 'none') +
#   geom_point() +
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#                      axis.title.x = element_text(size = 10),
#                      plot.title = element_text(face = "bold", size = 11)) +
#   scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) + 
#   ggtitle("MSSNG") +
#   xlab("Size Bin Start (bp)") + ylab("OR")
# 
# SSC.plot <- ggplot(SSC.less.30k, aes(x = size.bin.start, y = estimate)) +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high, alpha = 0.3))+ 
#   coord_cartesian(ylim=c(0,5)) +
#   scale_alpha(guide = 'none') + 
#   geom_point() +
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
#                      axis.title.x = element_text(size = 10),
#                      plot.title = element_text(face = "bold")) +
#   scale_x_continuous(breaks= c(0,5000,10000,15000,20000,25000,30000)) + 
#   ggtitle("SSC") +
#   xlab("Size Bin Start (bp)") + ylab("OR")
# 
# 
# ## Combine plots
# plot_grid(MSSNG.SSC.plot, MSSNG.plot, SSC.plot, nrow= 3)
# 
# plot.path <- "./figures/MSSNG.SSC.median.sizebins.png"
# ggsave(plot.path, width = 6, height = 7)

