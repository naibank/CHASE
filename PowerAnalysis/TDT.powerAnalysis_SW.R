################################################################################################################################################################################
# 
# TDT.powerAnalysis_SW.R
# purpose: plots power vs. sample size for TDT w/ current power labeled 
#         
# input: MSSNG.SSC_parent_proband_FisherTest_snvfreq0.01.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/
# output: MSSNG.SSC.TDT.missense.powerPlot.pRec0.9.snv0.01.png
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/powerAnalysis/figures/
#
# notes: power analysis only done on missense, pRec > 0.9, SNV freq < 0.01
#
##############################################################################################################################################################################

library(data.table)
library(reshape2)
library(ggplot2)
library(grid)
library(statmod)
library(pwr)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/powerAnalysis")

TDT.fisher.res <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data/MSSNG.SSC_parent_proband_FisherTest_snvfreq0.01.tsv", data.table = F)


Get_Power_Plot <- function(df){
  ## get power table for burden analysis or association analysis 
  pwr.n.table <- data.frame()
  
  df <- df[which(df$pRec == 0.9 & df$variant_type == "nonsynonymous"),]
  n <-sum(df$target.child, df$target.parent, df$bg.child, df$bg.parent) # total number of samples/individuals = 10437 
  n.var.per.sample <- n/10437
  
  # p.case = n1/(n1+n2)
  p.case <- df$target.child/(df$target.child + df$bg.child)
  
  # p.control = m2/(m1+m2 )
  p.control <- df$target.parent/(df$target.parent + df$bg.parent)
  
  h <- 2*asin(sqrt(p.case)) - 2*asin(sqrt(p.control))
  
  ### current power
  current.pwr <- (pwr.2p.test(h = h, n = n/2, sig.level = 0.05, alternative = "greater")$power) 
  current.pwr.row <- data.frame(power = current.pwr, N = n/n.var.per.sample, group = "current power")
  pwr.n.table <- rbind(pwr.n.table, current.pwr.row)
  
  ## power vs n data points
  for (pwr in seq(0.1, 0.9, 0.1)){
    n <-  pwr.2p.test(h = h, sig.level = 0.05, power = pwr, alternative = "greater")$n 
    n.for.pwr <- n*2/n.var.per.sample #round(n*2/3) ## number of families at __% power
    
    row <- data.frame(power = pwr, N = n.for.pwr, group = "calculated power")
    pwr.n.table <- rbind(pwr.n.table, row)
  }
  
  ## Plot power vs. n table
  plot <- ggplot() +
    geom_bar(data=subset(pwr.n.table, group == "calculated power"),
             aes(x = power, y = N), stat = "identity") + 
    geom_line(aes(x=c(current.pwr, current.pwr),y=c(0.2,10437)), linetype="dashed", colour = "red") + 
    geom_line(aes(x=c(current.pwr, 0),y=c(10437, 10437)), linetype="dashed", colour = "red") + 
    geom_point(data=subset(pwr.n.table, group== "current power"), 
               aes(x = power, y = N), colour = "red") +
    xlab("Power") + ylab("Number of Samples (n)") + 
    ggtitle("MSSNG+SSC TDT, Missense, pRec > 0.9, SNV_freq < 1%") +
    scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
    scale_y_continuous(breaks = seq(0, 60000, 10000)) +
    geom_text(aes(x=current.pwr, y=10437+6000, label='n = 10437 \n power = 0.44'), 
              color='red', size=5, fontface = "bold",) +
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15), 
          axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"))
  
  ggsave("./figures/MSSNG.SSC.TDT.missense.powerPlot.pRec0.9.snv0.01.png", plot = plot,
         width = 6, height = 5)
}

Get_Power_Plot(TDT.fisher.res)
