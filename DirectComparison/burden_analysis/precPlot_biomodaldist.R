setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

score <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

library(cutoff)
library('bbmle')
score <- score[!is.na(score$pRec), ]

mixmodel <- em(score$pRec, "normal","normal")
cutoff <- cutoff(mixmodel, distr = 2, type1 = 0) ##0.626

pdf("pRec.pdf", width = 4, height = 4)
hist(score$pRec, breaks = 100)
lines(x=c(cutoff[1],cutoff[1]), y=c(1,3000), lty = 2, col = "red")
text(round(cutoff[1], digits = 2), x=cutoff[2]-0.03, y = 1500, srt=90)
dev.off()
