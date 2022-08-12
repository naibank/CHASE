################################################################################################################################################################################
# 
# LogisticRegression_zscores_logtransformation_SW.R
#
# purpose: output logistic regression results on MSSNG+SSC z-score analysis 
# input: MSSNG_SSC_CH_hits_nofilt_sex_zscore.tsv
#           -> output of plot_logtransformation_zscores_SW.R
# output: e.g. logreg.zscores.pRec.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data
# notes:
#
##############################################################################################################################################################################

library(data.table)
library(dplyr)
library(broom)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/CH_count_del_size/data")

df.zscore <- as.data.frame(read.delim("./CH_hits/MSSNG_SSC_CH_hits_nofilt_sex_zscore.tsv"))
df.zscore$ASD.status <- ifelse(df.zscore$Relation %in% c("Proband-female", "Proband-male"), 1, 0)
df.zscore$sex <- ifelse(df.zscore$Relation %in% c("Proband-female", "Mother"), "female", "male")

df.zscore.pRec <- as.data.frame(read.delim("./CH_hits/MSSNG_SSC_CH_hits_pRec_sex_zscore.tsv"))
df.zscore.pRec$ASD.status <- ifelse(df.zscore.pRec$Relation %in% c("Proband-female", "Proband-male"), 1, 0)
df.zscore.pRec$sex <- ifelse(df.zscore.pRec$Relation %in% c("Proband-female", "Mother"), "female", "male")

Get_LogReg <- function(zscore.table, name){
  snv.count.res <- data.frame()
  for (SNV.count in c(0:8, "0-8", "1-8")){
    if (!SNV.count %in% c("0-8", "1-8")){
      df.count <- subset(zscore.table, zscore.table$CH_hit == SNV.count)
    }
    if (SNV.count == "0-8"){
      df.count <- zscore.table
    }
    if (SNV.count == "1-8"){
      df.count <- subset(zscore.table, CH_hit != 0)
    }
    formula <- ASD.status ~ sex + z.score 
    logit <- glm(formula, data = df.count, family = 'binomial')
    
    count <- data.frame(SNV.count = SNV.count)
    res.row <- broom::tidy(logit)[3,]
    res.row <- cbind(count, res.row)
    
    snv.count.res <- rbind(snv.count.res, res.row)
  }
  write.table(snv.count.res, sprintf("logreg.zscores.%s.tsv", name), 
              sep="\t", row.names=F, quote=F, col.names=T)
  return(snv.count.res)
}

logreg.nofilt <- Get_LogReg(df.zscore, "nofilt")
logreg.pRec <- Get_LogReg(df.zscore.pRec, "pRec")


# write.table(logreg.nofilt)


