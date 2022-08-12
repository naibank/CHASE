################################################################################################################################################################################
# 
# FisherTest_CNV_SV_rate_comparison_by_individual_SPARK_SW.R
# purpose: outputs SPARK Fisher's test results using spline of median, 
#           80th, & 90th percentiles for sliding window size bins
# input: SPARK_CH_hits (pRec > 0.9 and no pRec filt)
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/CH_hits
#        SPARK.2.5kb.n7.window.cutoffs.df.tsv (pRec > 0.9 and no pRec filt)
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/
# output: e.g. SPARK.fisher.spline.allcutoffs.tsv
#        /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/
# 
# notes: incomplete set of SPARK deletion data (only CNVs, not including SVs)
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

## Read CH hits and sliding window cut-offs
SPARK.CH.hits <- fread("./data/CH_hits/SPARK_CH_hits_nofilt_sex.tsv", data.table=F)
SPARK.CH.hits$Relation[which(SPARK.CH.hits$Relation %in% c("Proband-female", "Proband-male"))] <- "Proband"
SPARK.2.5kb.n7.window.cutoffs.df <- fread("./data/SPARK.2.5kb.n7.window.cutoffs.df.tsv", data.table=F)

SPARK.CH.hits.pRec <- fread("./data/CH_hits/SPARK_CH_hits_pRec_sex.tsv", data.table=F)
SPARK.CH.hits.pRec$Relation[which(SPARK.CH.hits.pRec$Relation %in% c("Proband-female", "Proband-male"))] <- "Proband"
SPARK.2.5kb.n7.pRec.window.cutoffs.df <- fread("./data/SPARK.2.5kb.n7.pRec.window.cutoffs.df", data.table=F)

########################################################################################################################################################################################
##### Median as Cut-off ####
Get_Spline_Median_Fisher_Test <- function(window.cutoffs.df, CH.hits.df){
  ### Return Fisher test results using control median slpine as cut-off.
  
  ## Take average of 20-30kb size bins as single bin
  window.cutoffs.df.20kb.less <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start <= 20000),]
  window.cutoffs.df.20kb.more <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start >= 20000),]
  window.cutoffs.df.20kb.less[41,]$size.bin.end <- 30000
  window.cutoffs.df.20kb.less[41,]$median <- mean(window.cutoffs.df.20kb.more$median)
  window.cutoffs.df.20kb.less[41,]$mean <- mean(window.cutoffs.df.20kb.more$mean)
  window.cutoffs.df.20kb.less[41,]$percentile_0.8 <- mean(window.cutoffs.df.20kb.more$percentile_0.8)
  window.cutoffs.df.20kb.less[41,]$percentile_0.9 <- mean(window.cutoffs.df.20kb.more$percentile_0.9)
  window.cutoffs.df.20kb.less[41,]$bin.count <- sum(window.cutoffs.df.20kb.more$bin.count)
  
  sliding.window.df <- window.cutoffs.df.20kb.less
  
  ## Remove > 30kb deletions and (0,0)
  n.30kb.more <- nrow(CH.hits.df[which(CH.hits.df$exSize > 30000),])
  CH.hits.df <- CH.hits.df[which(CH.hits.df$exSize <= 30000),]  # remove > 30kb
  CH.hits.df <- CH.hits.df[which(!(CH.hits.df$exSize == 0 & CH.hits.df$CH_hit == 0)),] # remove (0,0)
  
  ## Set del > 20kb to 20kb
  CH.hits.df[CH.hits.df$exSize > 20000, "exSize"]  <- 20000  
  
  ## Separate proband and parent CH hits
  CH.hits.df.probands <- CH.hits.df[which(CH.hits.df$Relation == "Proband"),]
  CH.hits.df.parents <- CH.hits.df[which(!CH.hits.df$Relation == "Proband"),]
  
  ## Get cut-offs from spline for Fisher test (expected # CH hits for each del size)
  spline.probands.cor <- spline(x = sliding.window.df$size.bin.start, 
                                y = sliding.window.df$median, 
                                xout = CH.hits.df.probands$exSize)
  CH.hits.df.probands$cutoff_median <- round(spline.probands.cor$y)
      # table(CH.hits.df.probands$CH_hit, CH.hits.df.probands$cutoff_median)
  
  spline.parents.cor <- spline(x = sliding.window.df$size.bin.start, 
                               y = sliding.window.df$median,
                               xout = CH.hits.df.parents$exSize)
  CH.hits.df.parents$cutoff_median <- round(spline.parents.cor$y)
  
  ## Find count of individuals above/below cut-off
  CH.hits.df.probands$cutoff_median_resi <- CH.hits.df.probands$CH_hit - CH.hits.df.probands$cutoff_median
  CH.hits.df.parents$cutoff_median_resi <- CH.hits.df.parents$CH_hit - CH.hits.df.parents$cutoff_median
  
  proband_above_expect <- sum(CH.hits.df.probands$cutoff_median_resi > 0)
  proband_below_expect <- sum(CH.hits.df.probands$cutoff_median_resi <= 0)
  control_above_expect <- sum(CH.hits.df.parents$cutoff_median_resi > 0)
  control_below_expect <- sum(CH.hits.df.parents$cutoff_median_resi <= 0)
  
  ## Fisher's test
  case_control_fisher_df <- data.frame(Above=c(proband_above_expect,control_above_expect),
                                       Below=c(proband_below_expect,control_below_expect))
  rownames(case_control_fisher_df) <- c('Case','Control')
  case_control_fisher_res <- fisher.test(case_control_fisher_df, alternative='greater')
  case_control_fisher_res.df <- broom::tidy(case_control_fisher_res)
  ## add counts to fisher results
  case_control_fisher_res.df$case.above <- proband_above_expect
  case_control_fisher_res.df$case.below <- proband_below_expect
  case_control_fisher_res.df$control.above <- control_above_expect
  case_control_fisher_res.df$control.below <- control_below_expect
  
  case_control_fisher_res.df[,1:3] <- signif(case_control_fisher_res.df[,1:3], 4) # round results
  
  ## add >30kb counts (removed from analysis)
  case_control_fisher_res.df$`n.>30kb (removed)`<- n.30kb.more
  
  return(case_control_fisher_res.df)
}

SPARK.median.fisher <- Get_Spline_Median_Fisher_Test(
  SPARK.2.5kb.n7.window.cutoffs.df, SPARK.CH.hits)
SPARK.median.fisher.10kb.more <- Get_Spline_Median_Fisher_Test(
  SPARK.2.5kb.n7.window.cutoffs.df[which(SPARK.2.5kb.n7.window.cutoffs.df$size.bin.start >= 10000),],
  SPARK.CH.hits[which(SPARK.CH.hits$exSize >= 10000),])

fisher.median <- rbind(SPARK.median.fisher, SPARK.median.fisher.10kb.more)
new.cols <- data.frame("cut-off" = c("median", "median"),
                       "deletion size range" = c("<30kb", "10-30kb"))
fisher.median <- cbind(new.cols, fisher.median[, c(1:4, 7:11)])
names(fisher.median)[names(fisher.median) == 'estimate'] <- 'OR'

## pRec > 0.9
SPARK.median.fisher.pRec <- Get_Spline_Median_Fisher_Test(
  SPARK.2.5kb.n7.pRec.window.cutoffs.df, SPARK.CH.hits.pRec )
SPARK.median.fisher.10kb.more.pRec <- Get_Spline_Median_Fisher_Test(
  SPARK.2.5kb.n7.pRec.window.cutoffs.df[which(SPARK.2.5kb.n7.pRec.window.cutoffs.df$size.bin.start >= 10000),],
  SPARK.CH.hits.pRec[which(SPARK.CH.hits.pRec$exSize >= 10000),])

fisher.median.pRec <- rbind(SPARK.median.fisher.pRec, SPARK.median.fisher.10kb.more.pRec)
new.cols <- data.frame("cut-off" = c("median", "median"),
                       "deletion size range" = c("<30kb", "10-30kb"))
fisher.median.pRec <- cbind(new.cols, fisher.median.pRec[, c(1:4, 7:11)])
names(fisher.median.pRec)[names(fisher.median.pRec) == 'estimate'] <- 'OR'


########################################################################################################################################################################################

##### 80th Percentile as Cut-off ####
Get_Spline_0.8_Fisher_Test <- function(window.cutoffs.df, CH.hits.df){
  ### Return Fisher test results using 80th percentile slpine as cut-off.
  
  ## Take average of 20-30kb size bins as single bin
  window.cutoffs.df.20kb.less <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start <= 20000),]
  window.cutoffs.df.20kb.more <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start >= 20000),]
  window.cutoffs.df.20kb.less[41,]$size.bin.end <- 30000
  window.cutoffs.df.20kb.less[41,]$median <- mean(window.cutoffs.df.20kb.more$median)
  window.cutoffs.df.20kb.less[41,]$mean <- mean(window.cutoffs.df.20kb.more$mean)
  window.cutoffs.df.20kb.less[41,]$percentile_0.8 <- mean(window.cutoffs.df.20kb.more$percentile_0.8)
  window.cutoffs.df.20kb.less[41,]$percentile_0.9 <- mean(window.cutoffs.df.20kb.more$percentile_0.9)
  window.cutoffs.df.20kb.less[41,]$bin.count <- sum(window.cutoffs.df.20kb.more$bin.count)
  
  sliding.window.df <- window.cutoffs.df.20kb.less
  
  ## Remove > 30kb deletions and (0,0)   
  n.30kb.more <- nrow(CH.hits.df[which(CH.hits.df$exSize > 30000),])
  CH.hits.df <- CH.hits.df[which(CH.hits.df$exSize <= 30000),]  # remove > 30kb
  CH.hits.df <- CH.hits.df[which(!(CH.hits.df$exSize == 0 & CH.hits.df$CH_hit == 0)),] # remove (0,0)
  
  ## Set del > 20kb to 20kb
  CH.hits.df[CH.hits.df$exSize > 20000, "exSize"]  <- 20000  
  
  ## Separate proband and parent CH hits
  CH.hits.df.probands <- CH.hits.df[which(CH.hits.df$Relation == "Proband"),]
  CH.hits.df.parents <- CH.hits.df[which(!CH.hits.df$Relation == "Proband"),]
  
  ## Get cut-offs from spline for Fisher test (expected # CH hits for each del size)
  spline.probands.cor <- spline(x = sliding.window.df$size.bin.start, 
                                y = sliding.window.df$percentile_0.8, 
                                xout = CH.hits.df.probands$exSize)
  CH.hits.df.probands$cutoff_percentile_0.8 <- round(spline.probands.cor$y)
  
  spline.parents.cor <- spline(x = sliding.window.df$size.bin.start, 
                               y = sliding.window.df$percentile_0.8,
                               xout = CH.hits.df.parents$exSize)
  CH.hits.df.parents$cutoff_percentile_0.8 <- round(spline.parents.cor$y)
  
  ## Find count of individuals above/below cut-off
  CH.hits.df.probands$cutoff_percentile_0.8_resi <- CH.hits.df.probands$CH_hit - CH.hits.df.probands$cutoff_percentile_0.8
  CH.hits.df.parents$cutoff_percentile_0.8_resi <- CH.hits.df.parents$CH_hit - CH.hits.df.parents$cutoff_percentile_0.8
  
  proband_above_expect <- sum(CH.hits.df.probands$cutoff_percentile_0.8_resi > 0)
  proband_below_expect <- sum(CH.hits.df.probands$cutoff_percentile_0.8_resi <= 0)
  control_above_expect <- sum(CH.hits.df.parents$cutoff_percentile_0.8_resi > 0)
  control_below_expect <- sum(CH.hits.df.parents$cutoff_percentile_0.8_resi <= 0)
  
  ## Fisher's test
  case_control_fisher_df <- data.frame(Above=c(proband_above_expect,control_above_expect),
                                       Below=c(proband_below_expect,control_below_expect))
  rownames(case_control_fisher_df) <- c('Case','Control')
  case_control_fisher_res <- fisher.test(case_control_fisher_df, alternative='greater')
  case_control_fisher_res.df <- broom::tidy(case_control_fisher_res)
  ## add counts to fisher results
  case_control_fisher_res.df$case.above <- proband_above_expect
  case_control_fisher_res.df$case.below <- proband_below_expect
  case_control_fisher_res.df$control.above <- control_above_expect
  case_control_fisher_res.df$control.below <- control_below_expect
  
  case_control_fisher_res.df[,1:3] <- signif(case_control_fisher_res.df[,1:3], 4) # round results
  
  ## add >30kb counts (removed from analysis)
  case_control_fisher_res.df$`n.>30kb (removed)`<- n.30kb.more
  
  return(case_control_fisher_res.df)
}

SPARK.percentile_0.8.fisher <- Get_Spline_0.8_Fisher_Test(
  SPARK.2.5kb.n7.window.cutoffs.df, SPARK.CH.hits)

SPARK.percentile_0.8.fisher.10kb.more <- Get_Spline_0.8_Fisher_Test(
  SPARK.2.5kb.n7.window.cutoffs.df[which(SPARK.2.5kb.n7.window.cutoffs.df$size.bin.start >= 10000),],
  SPARK.CH.hits[which(SPARK.CH.hits$exSize >= 10000),])

fisher.percentile_0.8 <- rbind(SPARK.percentile_0.8.fisher, SPARK.percentile_0.8.fisher.10kb.more)
new.cols <- data.frame("cut-off" = c("80th percentile", "80th percentile"),
                       "deletion size range" = c("<30kb", "10-30kb"))
fisher.percentile_0.8 <- cbind(new.cols, fisher.percentile_0.8[, c(1:4, 7:11)])
names(fisher.percentile_0.8)[names(fisher.percentile_0.8) == 'estimate'] <- 'OR'

## pRec > 0.9
SPARK.percentile_0.8.fisher.pRec <- Get_Spline_0.8_Fisher_Test(
  SPARK.2.5kb.n7.pRec.window.cutoffs.df, SPARK.CH.hits.pRec )
SPARK.percentile_0.8.fisher.10kb.more.pRec <- Get_Spline_0.8_Fisher_Test(
  SPARK.2.5kb.n7.pRec.window.cutoffs.df[which(SPARK.2.5kb.n7.pRec.window.cutoffs.df$size.bin.start >= 10000),],
  SPARK.CH.hits.pRec[which(SPARK.CH.hits.pRec$exSize >= 10000),])

fisher.percentile_0.8.pRec <- rbind(SPARK.percentile_0.8.fisher.pRec, SPARK.percentile_0.8.fisher.10kb.more.pRec)
new.cols <- data.frame("cut-off" = c("percentile_0.8", "percentile_0.8"),
                       "deletion size range" = c("<30kb", "10-30kb"))
fisher.percentile_0.8.pRec <- cbind(new.cols, fisher.percentile_0.8.pRec[, c(1:4, 7:11)])
names(fisher.percentile_0.8.pRec)[names(fisher.percentile_0.8.pRec) == 'estimate'] <- 'OR'


########################################################################################################################################################################################
##### 90th Percentile as Cut-off ####
Get_Spline_0.9_Fisher_Test <- function(window.cutoffs.df, CH.hits.df){
  ### Return Fisher test results using 90th percentile slpine as cut-off.
  
  ## Take average of 20-30kb size bins as single bin
  window.cutoffs.df.20kb.less <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start <= 20000),]
  window.cutoffs.df.20kb.more <- window.cutoffs.df[which(window.cutoffs.df$size.bin.start >= 20000),]
  window.cutoffs.df.20kb.less[41,]$size.bin.end <- 30000
  window.cutoffs.df.20kb.less[41,]$median <- mean(window.cutoffs.df.20kb.more$median)
  window.cutoffs.df.20kb.less[41,]$mean <- mean(window.cutoffs.df.20kb.more$mean)
  window.cutoffs.df.20kb.less[41,]$percentile_0.8 <- mean(window.cutoffs.df.20kb.more$percentile_0.8)
  window.cutoffs.df.20kb.less[41,]$percentile_0.9 <- mean(window.cutoffs.df.20kb.more$percentile_0.9)
  window.cutoffs.df.20kb.less[41,]$bin.count <- sum(window.cutoffs.df.20kb.more$bin.count)
  
  sliding.window.df <- window.cutoffs.df.20kb.less
  
  ## Remove > 30kb deletions and (0,0)   
  n.30kb.more <- nrow(CH.hits.df[which(CH.hits.df$exSize > 30000),])
  CH.hits.df <- CH.hits.df[which(CH.hits.df$exSize <= 30000),]  # remove > 30kb
  CH.hits.df <- CH.hits.df[which(!(CH.hits.df$exSize == 0 & CH.hits.df$CH_hit == 0)),] # remove (0,0)
  
  ## Set del > 20kb to 20kb
  CH.hits.df[CH.hits.df$exSize > 20000, "exSize"]  <- 20000  
  
  ## Separate proband and parent CH hits
  CH.hits.df.probands <- CH.hits.df[which(CH.hits.df$Relation == "Proband"),]
  CH.hits.df.parents <- CH.hits.df[which(!CH.hits.df$Relation == "Proband"),]
  
  ## Get cut-offs from spline for Fisher test (expected # CH hits for each del size)
  spline.probands.cor <- spline(x = sliding.window.df$size.bin.start, 
                                y = sliding.window.df$percentile_0.9, 
                                xout = CH.hits.df.probands$exSize)
  CH.hits.df.probands$cutoff_percentile_0.9 <- round(spline.probands.cor$y)
  
  spline.parents.cor <- spline(x = sliding.window.df$size.bin.start, 
                               y = sliding.window.df$percentile_0.9,
                               xout = CH.hits.df.parents$exSize)
  CH.hits.df.parents$cutoff_percentile_0.9 <- round(spline.parents.cor$y)
  
  ## Find count of individuals above/below cut-off
  CH.hits.df.probands$cutoff_percentile_0.9_resi <- CH.hits.df.probands$CH_hit - CH.hits.df.probands$cutoff_percentile_0.9
  CH.hits.df.parents$cutoff_percentile_0.9_resi <- CH.hits.df.parents$CH_hit - CH.hits.df.parents$cutoff_percentile_0.9
  
  proband_above_expect <- sum(CH.hits.df.probands$cutoff_percentile_0.9_resi > 0)
  proband_below_expect <- sum(CH.hits.df.probands$cutoff_percentile_0.9_resi <= 0)
  control_above_expect <- sum(CH.hits.df.parents$cutoff_percentile_0.9_resi > 0)
  control_below_expect <- sum(CH.hits.df.parents$cutoff_percentile_0.9_resi <= 0)
  
  ## Fisher's test
  case_control_fisher_df <- data.frame(Above=c(proband_above_expect,control_above_expect),
                                       Below=c(proband_below_expect,control_below_expect))
  rownames(case_control_fisher_df) <- c('Case','Control')
  case_control_fisher_res <- fisher.test(case_control_fisher_df, alternative='greater')
  case_control_fisher_res.df <- broom::tidy(case_control_fisher_res)
  ## add counts to fisher results
  case_control_fisher_res.df$case.above <- proband_above_expect
  case_control_fisher_res.df$case.below <- proband_below_expect
  case_control_fisher_res.df$control.above <- control_above_expect
  case_control_fisher_res.df$control.below <- control_below_expect
  
  case_control_fisher_res.df[,1:3] <- signif(case_control_fisher_res.df[,1:3], 4) # round results
  
  ## add >30kb counts (removed from analysis)
  case_control_fisher_res.df$`n.>30kb (removed)`<- n.30kb.more
  
  return(case_control_fisher_res.df)
}

SPARK.percentile_0.9.fisher <- Get_Spline_0.9_Fisher_Test(
  SPARK.2.5kb.n7.window.cutoffs.df, SPARK.CH.hits)

SPARK.percentile_0.9.fisher.10kb.more <- Get_Spline_0.9_Fisher_Test(
  SPARK.2.5kb.n7.window.cutoffs.df[which(SPARK.2.5kb.n7.window.cutoffs.df$size.bin.start >= 10000),],
  SPARK.CH.hits[which(SPARK.CH.hits$exSize >= 10000),])

fisher.percentile_0.9 <- rbind(SPARK.percentile_0.9.fisher, SPARK.percentile_0.9.fisher.10kb.more)
new.cols <- data.frame("cut-off" = c("90th percentile", "90th percentile"),
                       "deletion size range" = c("<30kb", "10-30kb"))
fisher.percentile_0.9 <- cbind(new.cols, fisher.percentile_0.9[, c(1:4, 7:11)])
names(fisher.percentile_0.9)[names(fisher.percentile_0.9) == 'estimate'] <- 'OR'

## pRec > 0.9
SPARK.percentile_0.9.fisher.pRec <- Get_Spline_0.9_Fisher_Test(
  SPARK.2.5kb.n7.pRec.window.cutoffs.df, SPARK.CH.hits.pRec )
SPARK.percentile_0.9.fisher.10kb.more.pRec <- Get_Spline_0.9_Fisher_Test(
  SPARK.2.5kb.n7.pRec.window.cutoffs.df[which(SPARK.2.5kb.n7.pRec.window.cutoffs.df$size.bin.start >= 10000),],
  SPARK.CH.hits.pRec[which(SPARK.CH.hits.pRec$exSize >= 10000),])

fisher.percentile_0.9.pRec <- rbind(SPARK.percentile_0.9.fisher.pRec, SPARK.percentile_0.9.fisher.10kb.more.pRec)
new.cols <- data.frame("cut-off" = c("percentile_0.9", "percentile_0.9"),
                       "deletion size range" = c("<30kb", "10-30kb"))
fisher.percentile_0.9.pRec <- cbind(new.cols, fisher.percentile_0.9.pRec[, c(1:4, 7:11)])
names(fisher.percentile_0.9.pRec)[names(fisher.percentile_0.9.pRec) == 'estimate'] <- 'OR'



df.out <- rbind(fisher.median, fisher.percentile_0.8, fisher.percentile_0.9)
write.table(df.out, "./data/SPARK.fisher.spline.allcutoffs.tsv", 
            sep="\t", row.names=F, quote=F, col.names=T)

df.out.pRec <- rbind(fisher.median.pRec, fisher.percentile_0.8.pRec, fisher.percentile_0.9.pRec)
write.table(df.out.pRec, "./data/SPARK.fisher.spline.allcutoffs.pRec.tsv", 
            sep="\t", row.names=F, quote=F, col.names=T)
