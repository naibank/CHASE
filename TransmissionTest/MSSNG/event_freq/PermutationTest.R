set.seed(13)
library(data.table)
library(ggplot2)

perm_score <- function(scores, len){
  rand_scores <- c()
  for(i in 1:10000){
    rand_scores <- c(rand_scores, sum(sample(scores, len), na.rm = T))
  }
  
  return(rand_scores)
}

genemania <- read.delim("../../MSSNG/network_analysis/MSSNG_sanitized.txt-results.scores.txt",
                        header = F, col.names = c("gene", "score"))

genemania <- genemania[!is.na(genemania$score), ]

genes <- data.table::fread("../../../../../ReferenceData/geneInfo2019/hg38_refGene_20200708.transcript.txt", data.table = F)
genes <- unique(genes[, c("V5", "V6")])

genemania <- merge(genemania, genes, by.x = "gene", by.y = "V5", all = F)
names(genemania)[3] <- c("enzid")
# gsMain <- gsMain[1:36]
load("../../gene-sets/gsMain_PGC_2021.RData")

genemania <- genemania[genemania$enzid %in% unlist(gsMain), ]
gsMain <- append(gsMain[18:24], gsMain[31:32]) # extract control gene sets


test.out <- data.table()

for (cutoff in c(seq(0, 0.1, by = 0.01))){
  genemania <- genemania[genemania$score > cutoff, ]
  
  for (i in 1:length(gsMain)){
    scoresum <- sum(genemania$score[genemania$enzid %in% gsMain[[i]]], na.rm = T)
    len <- sum(genemania$enzid %in% gsMain[[i]])
    ctrl <- perm_score(genemania$score, len)
    
    # p <- t.test(x=ctrl, mu=scoresum, alternative = "less")
    # p$p.value
    
    p <- sum(ctrl >= scoresum)/length(ctrl)
    
    test.out <- rbind(test.out, data.frame("score_cut-off" = cutoff,
                                           "gene-set" = names(gsMain[i]),
                                           "P" = p))
  }
}

write.table(test.out, "permutation_cutoff_optimization.tsv", row.names=F, quote=F, col.names=T)


## plot false discovery rate
FPR <- data.frame("score_cut-off" = c(seq(0, 0.1, by = 0.01)),
                  "FPR" = NA)

for (i in 1:nrow(FPR)){
  cutoff <- FPR[i,]$score_cut.off
  tmp <- test.out[which(test.out$score_cut.off == cutoff),]
  nFP <- nrow(tmp[which(tmp$P < 0.05),])
  n_geneset <- 9
  FPR[i,]$FPR <- signif(nFP/n_geneset, 2)
}

ggplot(FPR, aes(x = factor(score_cut.off), y = FPR)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_label(aes(label = FPR)) + 
  xlab("Score Cut-off") + ylab("False Positive Rate") 

ggsave("permutation_cutoff_optimization.png", height = 5, width = 7)


###### test MSSNG
load("../../gene-sets/gsMain_PGC_2021.RData")

genemania <- read.delim("../../MSSNG/network_analysis/MSSNG_sanitized.txt-results.scores.txt",
                        header = F, col.names = c("gene", "score"))

genemania <- genemania[!is.na(genemania$score), ]

genemania <- merge(genemania, genes, by.x = "gene", by.y = "V5", all = F)
names(genemania)[3] <- c("enzid")
genemania <- genemania[genemania$enzid %in% unlist(gsMain), ]

gsMain <- gsMain[1:36]

genemania <- genemania[genemania$score > 0.01, ]
test.out <- data.table()

for (i in 1:length(gsMain)){
  scoresum <- sum(genemania$score[genemania$enzid %in% gsMain[[i]]], na.rm = T)
  len <- sum(genemania$enzid %in% gsMain[[i]])
  ctrl <- perm_score(genemania$score, len)
  
  # p <- t.test(x=ctrl, mu=scoresum, alternative = "less")
  # p$p.value
  
  p <- sum(ctrl >= scoresum)/length(ctrl)
  
  test.out <- rbind(test.out, data.frame("gene-set" = names(gsMain[i]),
                                         "P" = p))
}
mssng <- test.out[order(test.out$P), ]

###### test SSC unaff
load("../../gene-sets/gsMain_PGC_2021.RData")

genemania <- read.delim("../../SSC/network_analysis/SSC_unaffectedsiblings_sanitized.txt-results.scores.txt",
                        header = F, col.names = c("gene", "score"))

genemania <- genemania[!is.na(genemania$score), ]

genemania <- merge(genemania, genes, by.x = "gene", by.y = "V5", all = F)
names(genemania)[3] <- c("enzid")
genemania <- genemania[genemania$enzid %in% unlist(gsMain), ]

gsMain <- gsMain[1:36]

genemania <- genemania[genemania$score > 0.01, ]
test.out <- data.table()

for (i in 1:length(gsMain)){
  scoresum <- sum(genemania$score[genemania$enzid %in% gsMain[[i]]], na.rm = T)
  len <- sum(genemania$enzid %in% gsMain[[i]])
  ctrl <- perm_score(genemania$score, len)
  
  # p <- t.test(x=ctrl, mu=scoresum, alternative = "less")
  # p$p.value
  
  p <- sum(ctrl >= scoresum)/length(ctrl)
  
  test.out <- rbind(test.out, data.frame("gene-set" = names(gsMain[i]),
                                         "P" = p))
}

unaff.ssc <- test.out[order(test.out$P), ]


###### test SSC proband
load("../../gene-sets/gsMain_PGC_2021.RData")

genemania <- read.delim("../../SSC/network_analysis/SSC_probands_sanitized.txt-results.scores.txt",
                        header = F, col.names = c("gene", "score"))

genemania <- genemania[!is.na(genemania$score), ]

genemania <- merge(genemania, genes, by.x = "gene", by.y = "V5", all = F)
names(genemania)[3] <- c("enzid")
genemania <- genemania[genemania$enzid %in% unlist(gsMain), ]

gsMain <- gsMain[1:36]

genemania <- genemania[genemania$score > 0.01, ]
test.out <- data.table()

for (i in 1:length(gsMain)){
  scoresum <- sum(genemania$score[genemania$enzid %in% gsMain[[i]]], na.rm = T)
  len <- sum(genemania$enzid %in% gsMain[[i]])
  ctrl <- perm_score(genemania$score, len)
  
  # p <- t.test(x=ctrl, mu=scoresum, alternative = "less")
  # p$p.value
  
  p <- sum(ctrl >= scoresum)/length(ctrl)
  
  test.out <- rbind(test.out, data.frame("gene-set" = names(gsMain[i]),
                                         "P" = p))
}

prob.ssc <- test.out[order(test.out$P), ]

mssng$set <- "MSSNG affected"
unaff.ssc$set <- "SSC unaffected"
prob.ssc$set <- "SSC affected"

dt.plot <- rbind(mssng, rbind(prob.ssc, unaff.ssc))
dt.plot$BHFDR <- NA
for(set in unique(dt.plot$set)){
  dt.plot$BHFDR[dt.plot$set == set] <- p.adjust(dt.plot$P[dt.plot$set == set], method = "BH")
}

dt.plot$gene.set <- factor(dt.plot$gene.set, levels = rev(unique(dt.plot$gene.set[order(dt.plot$BHFDR)])))
library(ggplot2)
ggplot(dt.plot, aes(y = gene.set, x = -log10(BHFDR))) +
  geom_bar(stat = "identity") + facet_wrap(set~., nrow = 1) + theme_bw() +
  geom_vline(xintercept = -log10(0.1), lty = 2, color = "red")

ggsave("gene.mania.by.score.png", width = 11, height = 7)
