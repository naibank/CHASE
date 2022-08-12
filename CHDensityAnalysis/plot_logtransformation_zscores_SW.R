################################################################################################################################################################################
# 
# plot_logtransformation_zscores_SW.R
# purpose: - outputs MSSNG+SSC CH hits table with zscore 
#          - plots SNV counts vs. Z-score plot for MSSNG+SSC 
# input: MSSNG+SSC CH hits tables with sex: MSSNG_SSC_CH_hits_nofilt_sex.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data
# output: - CH hits table with zscore: MSSNG_SSC_CH_hits_pRec_sex_zscore.tsv
#               -> used as input for LogisticRegression_zscores_logtransformation_SW.R 
#         - SNV counts vs. Z-score plot: MSSNG.SSC.zscores.nofilt.png
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/figures
# notes:
#
##############################################################################################################################################################################


library(data.table)
library(dplyr)
library(ggplot2)
library(vioplot)
library(forcats)
library(hrbrthemes)
library(viridis)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette.t <- c("#99999933", "#E69F0033", "#56B4E933", "#009E7333", "#F0E44233", "#0072B233", "#D55E0033", "#CC79A733")
cbVermeer <- c('#6495ED','#93CCEA','#6495ED50','#93CCEA50','#6495ED20','#93CCEA20')

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data")

df.ch <- as.data.frame(read.delim("./CH_hits/MSSNG_SSC_CH_hits_nofilt_sex.tsv"))
df.ch <- subset(df.ch,df.ch$exSize>0)

df.ch.pRec <- as.data.frame(read.delim("./CH_hits/MSSNG_SSC_CH_hits_pRec_sex.tsv"))
df.ch.pRec <- subset(df.ch.pRec,df.ch.pRec$exSize>0)


## Plot without log transformation
plot(log(df.ch$exSize),df.ch$CH_hit)
plot(df.ch$CH_hit, df.ch$exSize, ylim=c(0,50000),xlim=c(0,15),
     ylab="length",xlab="SNV counts")
for (i in 0:10) {
  vioplot(subset(df.ch$exSize,df.ch$CH_hit ==i), at = i, col=cbPalette.t[6],wex=.5, add=T)
}
median(subset(df.ch$exSize,df.ch$CH_hit ==8))
lines(c(0,8), c(782,11100),col=cbPalette[8],lwd=4,lty=2)

## Plot log transformed length vs. count
plot(log(df.ch$exSize),log(df.ch$CH_hit+1), axes=F, xlab="log length",ylab ='log count+1',
     pch=19,cex=.76,cex.lab=1.5)
axis(side=1,at=log(c(1,seq(25,250,25),seq(500,1000,250),seq(2000,20000,2000))),c(0,seq(25,250,25),seq(500,1000,250),seq(2000,20000,2000)))
axis(side=2,at=log(seq(1,15,1)),seq(1,15,1))
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==0)), horizontal = T, at = 0, add=T,wex=.25, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==1)), horizontal = T, at = log(2), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==2)), horizontal = T, at = log(3), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==3)), horizontal = T, at = log(4), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==4)), horizontal = T, at = log(5), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==5)), horizontal = T, at = log(6), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==6)), horizontal = T, at = log(7), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==7)), horizontal = T, at = log(8), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==8)), horizontal = T, at = log(9), add=T,wex=.15, col=cbPalette.t[6])

## Plot log transformed count vs. length
plot(log(df.ch$CH_hit+1),log(df.ch$exSize), axes=F, ylab="log length",xlab ='log count+1',
     pch=19,cex=.76,cex.lab=1.5)
axis(side=2,at=log(c(1,seq(25,250,25),seq(500,1000,250),seq(2000,20000,2000))),c(0,seq(25,250,25),seq(500,1000,250),seq(2000,20000,2000)))
axis(side=1,at=log(seq(1,15,1)),seq(1,15,1))
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==0)), horizontal = F, at = 0, add=T,wex=.25, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==1)), horizontal = F, at = log(2), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==2)), horizontal = F, at = log(3), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==3)), horizontal = F, at = log(4), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==4)), horizontal = F, at = log(5), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==5)), horizontal = F, at = log(6), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==6)), horizontal = F, at = log(7), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==7)), horizontal = F, at = log(8), add=T,wex=.15, col=cbPalette.t[6])
vioplot(log(subset(df.ch$exSize,df.ch$CH_hit ==8)), horizontal = F, at = log(9), add=T,wex=.15, col=cbPalette.t[6])

lines(c(0,log(15)),c(log(750),log(20000)),lwd=3,lty=2,col=cbPalette[7])



#### Plot SNV count vs. Z-score (of log length) ####

## Add col for z-scores
Get_Zscore <- function(CH_hits_table){
  df.zscore <- data.frame()
  for (SNV.count in 0:8){
    if (SNV.count == 8){
      CH_hits_table.count <- subset(CH_hits_table, CH_hits_table$CH_hit > 7) # 8+ SNV count
      CH_hits_table.count$CH_hit <- 8 # mutate CH_hit > 7 to 8
    } else {
      CH_hits_table.count <- subset(CH_hits_table, CH_hits_table$CH_hit == SNV.count)
    }
    CH_hits_table.count$log.exSize <- log(CH_hits_table.count$exSize)
    CH_hits_table.count$z.score <- ((CH_hits_table.count$log.exSize - mean(CH_hits_table.count$log.exSize))/
                              sd(CH_hits_table.count$log.exSize))
    
    df.zscore <- rbind(df.zscore, CH_hits_table.count)
  }
  return(df.zscore)
}

zscore.table.nofilt <- Get_Zscore(df.ch)
zscore.table.pRec <- Get_Zscore(df.ch.pRec)

# write.table(zscore.table.nofilt, "./CH_hits/MSSNG_SSC_CH_hits_nofilt_sex_zscore.tsv", 
#             sep="\t", row.names=F, quote=F, col.names=T)
write.table(zscore.table.pRec, "./CH_hits/MSSNG_SSC_CH_hits_pRec_sex_zscore.tsv",
            sep="\t", row.names=F, quote=F, col.names=T)

## Plot Z-scores 
zscore.table.nofilt$CH_hit <- as.factor(zscore.table.nofilt$CH_hit)

no.zero <- subset(zscore.table.nofilt, CH_hit != 0)
no.zero$CH_hit <- 10

ggplot(zscore.table.nofilt, aes(y = z.score, x = CH_hit, fill = Relation)) + #zscore.table.nofilt,
  geom_violin(data =zscore.table.nofilt, position=position_dodge(width = .6),size= .4) +
  geom_boxplot(data=zscore.table.nofilt[-9183,], width=.2, alpha=0.5, position=position_dodge(width = .6),size= .4) +
  geom_violin(data= no.zero, position=position_dodge(width = .9),size= .4) +
  geom_boxplot(data=no.zero, width=.2, alpha=0.5, position=position_dodge(width = .9),size= .4) +
  scale_fill_brewer(palette = "Accent") +
  geom_hline(yintercept=0, linetype='dashed', col = 'red') +
  xlab("SNV count") + ylab("z-score") +
  scale_x_discrete(labels = c(" ", 0:7, "8+", "1-8"), limits = 0:10) 
  
ggsave("../figures/MSSNG.SSC.zscores.nofilt.png", width = 10, height = 3)

