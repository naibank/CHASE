################################################################################################################################################################################
# 
# 20220712_logtransformation_CHASE_Elemi.R
# purpose: Elemi's script for plotting log transformations for CH event density analysis
# input: MSSNG_SSC_CH_hits_nofilt.tsv
# output: (1) log-log plot, for both the deletion size and the SNV count
#         (2) same plot, but switching the y and x-axes
#         (3) plot without transformations
#         (4) z-scores for all SNVs counts, and plotted them a consecutive order for increasing SNVs
# 
# notes:
#
##############################################################################################################################################################################

setwd("C:/Users/elemi breetvelt/OneDrive - SickKids/Research/Projects/CHASE")

df.ch <- as.data.frame(read.delim("MSSNG_SSC_CH_hits_nofilt.tsv"))
df.ch <- subset(df.ch,df.ch$exSize>0)


plot(log(df.ch$exSize),df.ch$CH_hit)
plot(df.ch$CH_hit, df.ch$exSize, ylim=c(0,50000),xlim=c(0,15),
     ylab="length",xlab="SNV counts")
for (i in 0:10) {
  vioplot(subset(df.ch$exSize,df.ch$CH_hit ==i), at = i, col=cbPalette.t[6],wex=.5, add=T)
}
median(subset(df.ch$exSize,df.ch$CH_hit ==8))
lines(c(0,8), c(782,11100),col=cbPalette[8],lwd=4,lty=2)



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

yi <- log(subset(df.ch$exSize,df.ch$CH_hit ==0))
yin <- ((yi-mean(yi))/sd(yi))


vioplot(yin, at =1, horizontal = F,col=cbPalette.t[7],xlim=c(0,11),ylim=c(-4,4),
        axes=F,cex.axis=1.25)
axis(1,at = seq(1,9,1),c(seq(0,7,1),">=8"),cex.axis=1.25)
for (i in 1:7){
  yi <- log(subset(df.ch$exSize,df.ch$CH_hit ==i))
  yin <- ((yi-mean(yi))/sd(yi))  
  vioplot(yin, at =1+i, horizontal = F,col=cbPalette.t[7],add=T,axes=F)
}
yi <- log(subset(df.ch$exSize,df.ch$CH_hit > 7))
yin <- ((yi-mean(yi))/sd(yi))  
vioplot(yin, at =9, horizontal = F,col=cbPalette.t[7],add=T,axes=F)  

abline(h=0,lty=2,lwd=3,col='gray75')
text(.1,3,"z-score",srt=90,font=3,cex=1.25)
text(9,-3.85,"SNV count",srt=0,font=3,cex=1.25)




