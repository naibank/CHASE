################################################################################################################################################################################
# 
# BA.powerAnalysis_SW.R
# purpose: plots power vs. sample size for BA w/ current power labeled 
#         
# input: MSSNG.SSC_parent_proband.SNV1p.pRec0.9_clogit_res.tsv
#        /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_BA/data/LogRegResultsWithOR/
# output:MSSNG.SSC.BA.missense.powerPlot.pRec0.9.snv0.01.png
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/powerAnalysis/figures/
#
# notes: power analysis only done on missense, pRec > 0.9, SNV freq < 0.01
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

## read clogit (burden analysis) results
BA.clogit.res <- fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_BA/data/LogRegResultsWithOR/MSSNG.SSC_parent_proband.SNV1p.pRec0.9_clogit_res.tsv", data.table=F)
# MSSNG+SSC burden analysis, pRec > 0.9, no SNV freq filter; from unbiased_burden_analysis

Get_Power_Plot <- function(df, method = "BA"){
  ## get power table for burden analysis or association analysis 
  if (method == "BA"){
    pwr.n.table <- data.frame()
    missense.res <- df[which(df$variant.type == "mis"),]
    n <- 10437
    u <- 1 # only one coefficient - CH Freq
    w <- 2 # number of covariates - Sex & totalCHevents
    
    pR.ref = 1 - missense.res$ref.deviance / missense.res$ref.null.deviance # variance accounted for in the pop. by covariate sex
    pR.add = 1 - missense.res$add.deviance / missense.res$add.null.deviance # variance accounted for in the pop. by covariate sex & CHrate
    
    f2 <- (pR.add - pR.ref) / (1 - pR.add) # based on https://www.statmethods.net/stats/power.html two models
    
    ### current power
    current.pwr <- pwr.f2.test(u = u, v = n - w - u - 1, f2 = f2, sig.level = 0.05)$power
    current.pwr.row <- data.frame(power = current.pwr, N = n, group = "current power")
    pwr.n.table <- rbind(pwr.n.table, current.pwr.row)
    
    ## power vs n data points
    for (pwr in seq(0.1, 0.9, 0.1)){
      v <- pwr.f2.test(u = u, power = pwr, f2 = f2, sig.level = 0.05)$v # get sample size with pwr power
      n.for.pwr <- v + w + u + 1 # get n to reach specific power; v = n - w - u - 1
      
      row <- data.frame(power = pwr, N = n.for.pwr, group = "calculated power")
      pwr.n.table <- rbind(pwr.n.table, row)
    }
    
    ## Plot power vs. n table
    plot <- ggplot() +
      geom_bar(data=subset(pwr.n.table, group == "calculated power"),
               aes(x = power, y = N), stat = "identity") + 
      geom_line(aes(x=c(current.pwr, current.pwr),y=c(0.2,n)), linetype="dashed", colour = "red") + 
      geom_line(aes(x=c(current.pwr, 0),y=c(n,n)), linetype="dashed", colour = "red") + 
      geom_point(data=subset(pwr.n.table, group== "current power"), 
                 aes(x = power, y = N), colour = "red") +
      xlab("Power") + ylab("Number of Samples (n)") + 
      ggtitle("MSSNG+SSC Burden Analysis - Missense, pRec > 0.9, SNV freq < 1%") +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(breaks = seq(0, max(pwr.n.table$N) +50000, 50000))+
      geom_text(aes(x=current.pwr, y=n+45000,  fontface = "bold",
                    label= sprintf('n = 10437 \n power = %s', signif(current.pwr, 2))), 
                color='red', size=5) +
      theme(axis.text.x = element_text(size = 15), 
            axis.text.y = element_text(size = 15), 
            axis.title.y = element_text(size = 15, face = "bold"),
            axis.title.x = element_text(size = 15, face = "bold"))
    
    ggsave("./figures/MSSNG.SSC.BA.missense.powerPlot.pRec0.9.snv0.01.png", plot = plot,
           width = 6, height = 5)
  }
}

Get_Power_Plot(BA.clogit.res)


# ## example:
# dt <- read.delim("table.for.power.exonic.tsv", stringsAsFactors = F)
# n <- 474
# u <- 1
# w <- 2
# ref <- "status ~ sex + ncloci"
# add <- paste(ref, "+", "ExpansionCount")
# 
# lm.ref <- glm(ref, dt, family = binomial(link = "logit"))
# lm.add <- glm(add, dt, family = binomial(link = "logit"))
# 
# ano <- anova(lm.ref, lm.add, test = "Chisq")
# pvalue <- ano$`Pr(>Chi)`[2]
# s <- summary(lm.add)
# 
# pR1 = 1 - lm.ref$deviance / lm.ref$null.deviance
# pR2 = 1 - lm.add$deviance / lm.add$null.deviance # R2 based on deviance
# 
# f2 <- (pR2 - pR1) / (1 - pR2)#based on https://www.statmethods.net/stats/power.html two models
# 
# ###current power
# pwr.f2.test(u = u, v = n - w - u - 1, f2 = f2, sig.level = 0.05)$power
# 
# ### get sample size with 80% power
# pwr.f2.test(u = u, power = 0.8, f2 = f2, sig.level = 0.05)$v
# #v = n - w - u - 1
# n = 648.2805 + w + u + 1
# 
# ###########################################
# # obtain sample sizes
# # range of correlations
# r <- seq(0.005, 0.03, 0.0001)
# n <- c(474)
# # obtain sample sizes
# dt.tmp <- data.frame()
# for (i in n){
#   for (j in r){
#     result <- pwr.f2.test(u = u, v = i - u - w - 1, f2 = j,
#                          sig.level = .1)
#     dt.tmp <- rbind(dt.tmp, data.frame("Test" = "two-sided", "effectSize" = j, "power" = result$power, stringsAsFactors = F))
#     
#     # result <- pwr.f2.test(u = u, v = i - u - w - 1, f2 = j,
#     #                       sig.level = .05)
#     # dt.tmp <- rbind(dt.tmp, data.frame("Test" = "one-sided", "effectSize" = j, "power" = result$power, stringsAsFactors = F))
#   }
# }
# 
# # dt.tmp$sampleSize <- factor(dt.tmp$sampleSize)
# ggplot(dt.tmp, aes(x = effectSize, y = power, color = Test, group = Test)) + geom_line() + theme_bw() +
#   ggtitle("Power calculation at observed effect sizes (n = 474)") + geom_vline(xintercept = c(0.012, 0.02), lty = 2, lwd = .3)
# 
# ggsave("power.analysis.png", width = 7, height = 5)
