################################################################################################################################################################################
# 
# association.powerAnalysis_SW.R
# purpose: plots power vs. sample size for spline CH density analysis w/ current power labeled 
#         
# input: MSSNG.SSC.fisher.spline.allcutoffs.pRec.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/
# output: MSSNG.SSC.association.90perc.10kb.powerPlot.pRec0.9.png
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/powerAnalysis/figures/
#
# notes: power analysis only done on 90th percentile, pRec > 0.9, del size > 10kb
#
##############################################################################################################################################################################


library(data.table)
library(reshape2)
library(ggplot2)
library(grid)
library(pwr)

# The numerator degrees of freedom (u) is number of coefficients in the full model excluding the intercept and covariates (as you mentioned earlier). 
# The denominator degrees of freedom (v) is equal to N - u - w - 1, where N is number of samples, and u is the numerator degrees of freedom, w is the number of covariates.


setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/powerAnalysis")

## read fisher (burden analysis) results
association.fisher.res <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data/MSSNG.SSC.fisher.spline.allcutoffs.pRec.tsv", data.table=F)
# MSSNG+SSC burden analysis, pRec > 0.9, no SNV freq filter; from unbiased_burden_analysis

Get_Power_Plot <- function(df){
  ## get power table for burden analysis or association analysis 
 
  pwr.n.table <- data.frame()
  target.res <- df[which(df$cut.off == "percentile_0.9" & df$deletion.size.range == "10-30kb"),] # 90th percentile, 10kb+
  
  n <-sum(target.res$case.above, target.res$case.below, target.res$control.above, target.res$control.below) # total number of samples/individuals = 10437 
  n.var.per.sample <- n/10437
  
  # p.case = n1/(n1+n2)
  p.case <- target.res$case.above/(target.res$case.above + target.res$case.below)
  
  # p.control = m2/(m1+m2 )
  p.control <- target.res$control.above/(target.res$control.above + target.res$control.below)
  
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
    ggtitle("MSSNG+SSC Spline Assoc., 90th percentile, 10kb+, pRec > 0.9") +
    scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
    scale_y_continuous(breaks = seq(0, max(pwr.n.table$N) +50000, 50000))+
    geom_text(aes(x=current.pwr, y=n+48000, 
                  label=sprintf('n = 10437 \n power = %s', signif(current.pwr, 2))), 
              color='red', size=5, fontface = "bold",) +
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15), 
          axis.title.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"))
    
    ggsave("./figures/MSSNG.SSC.association.90perc.10kb.powerPlot.pRec0.9.png", plot = plot,
           width = 6, height = 5)
}

Get_Power_Plot(association.fisher.res)
