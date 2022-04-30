library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

#Load("../Unbiased/unbiased_synchevent_covariate.RData")

Get_ConservationScore_Groups <- function(meta,
                               SNV_data,
                               proband_IDs,
                               father_IDs,
                               mother_IDs,
                               type='synonymous SNV',
                               score_thresh=0.1) {
  sj_df <- SNV_data[SNV_data$X.Sample %in% meta$Sample.ID
                    & SNV_data$effect_priority == type
                    & SNV_data$gnomAD_pRec >= 0.9 & !is.na(SNV_data$gnomAD_pRec)
                    & SNV_data$gnomAD_oe_lof_upper >= 0.35 & !is.na(SNV_data$gnomAD_oe_lof_upper)
                    ,]
  
  passed <- sj_df$X.Sample[sj_df$phylopMam_avg >= 1]
  sj_df$Pass <- FALSE
  sj_df$Pass[sj_df$X.Sample %in% passed] <- TRUE
  
  sj_df$Relation <- 'NA'
  sj_df$Relation[sj_df$X.Sample %in% proband_IDs] <- 'Proband'
  sj_df$Relation[sj_df$X.Sample %in% father_IDs] <- 'Father'
  sj_df$Relation[sj_df$X.Sample %in% mother_IDs] <- 'Mother'
  
  sj_df$Group <- 'NA'
  sj_df$Group[sj_df$X.Sample %in% proband_IDs] <- 'Child'
  sj_df$Group[sj_df$X.Sample %in% father_IDs |
                sj_df$X.Sample %in% mother_IDs] <- 'Parent'
  
  return (sj_df[,c('X.Sample','Relation', 'Group', 'Pass')])
}

res <- data.frame()

MSSNG_parent_proband_ConservationScore_Groups <- Get_ConservationScore_Groups(MSSNG_meta,
                                                          MSSNG_parent_proband_proc_SNVs,
                                                          MSSNG_proband_IDs,
                                                          MSSNG_father_IDs,
                                                          MSSNG_mother_IDs)
MSSNG_parent_proband_ConservationScore_Counts <- data.frame(MSSNG_parent_proband_ConservationScore_Groups %>% group_by(Group) %>%
                                                    summarise(Passed=sum(Pass), Failed=sum(!Pass)))
rownames(MSSNG_parent_proband_ConservationScore_Counts) <- MSSNG_parent_proband_ConservationScore_Counts$Group
MSSNG_parent_proband_ConservationScore_Counts <- MSSNG_parent_proband_ConservationScore_Counts[,-1]

#Fisher's exact test
MSSNG_parent_proband_ConservationScore_Fisher <- fisher.test(MSSNG_parent_proband_ConservationScore_Counts,
                                                   alternative='greater')
res <- rbind(res, data.frame(Name='MSSNG_parent_proband',
                             Child_Passed=MSSNG_parent_proband_ConservationScore_Counts$Passed[1],
                             Child_Failed=MSSNG_parent_proband_ConservationScore_Counts$Failed[1],
                             Parent_Passed=MSSNG_parent_proband_ConservationScore_Counts$Passed[2],
                             Parent_Failed=MSSNG_parent_proband_ConservationScore_Counts$Failed[2],
                             p_val=MSSNG_parent_proband_ConservationScore_Fisher$p.value,
                             OR=MSSNG_parent_proband_ConservationScore_Fisher$estimate))

MSSNG_parent_proband_0.1CH_ConservationScore_Groups <- Get_ConservationScore_Groups(MSSNG_meta,
                                                                MSSNG_parent_proband_0.1CH_SNVs,
                                                                MSSNG_proband_IDs,
                                                                MSSNG_father_IDs,
                                                                MSSNG_mother_IDs)

MSSNG_parent_proband_0.1CH_ConservationScore_Counts <- data.frame(MSSNG_parent_proband_0.1CH_ConservationScore_Groups %>% group_by(Group) %>%
                                                          summarise(Passed=sum(Pass), Failed=sum(!Pass)))
rownames(MSSNG_parent_proband_0.1CH_ConservationScore_Counts) <- MSSNG_parent_proband_0.1CH_ConservationScore_Counts$Group
MSSNG_parent_proband_0.1CH_ConservationScore_Counts <- MSSNG_parent_proband_0.1CH_ConservationScore_Counts[,-1]

#Fisher's exact test
MSSNG_parent_proband_0.1CH_ConservationScore_Fisher <- fisher.test(MSSNG_parent_proband_0.1CH_ConservationScore_Counts,
                                                         alternative='greater')

res <- rbind(res, data.frame(Name='MSSNG_parent_proband_0.1CH',
                             Child_Passed=MSSNG_parent_proband_0.1CH_ConservationScore_Counts$Passed[1],
                             Child_Failed=MSSNG_parent_proband_0.1CH_ConservationScore_Counts$Failed[1],
                             Parent_Passed=MSSNG_parent_proband_0.1CH_ConservationScore_Counts$Passed[2],
                             Parent_Failed=MSSNG_parent_proband_0.1CH_ConservationScore_Counts$Failed[2],
                             p_val=MSSNG_parent_proband_0.1CH_ConservationScore_Fisher$p.value,
                             OR=MSSNG_parent_proband_0.1CH_ConservationScore_Fisher$estimate))

SSC_parent_proband_0.01CH_ConservationScore_Groups <- Get_ConservationScore_Groups(SSC_meta,
                                                               SSC_parent_proband_0.01CH_SNVs,
                                                               SSC_proband_IDs,
                                                               SSC_father_IDs,
                                                               SSC_mother_IDs)

SSC_parent_proband_0.01CH_ConservationScore_Counts <- data.frame(SSC_parent_proband_0.01CH_ConservationScore_Groups %>% group_by(Group) %>%
                                                         summarise(Passed=sum(Pass), Failed=sum(!Pass)))
rownames(SSC_parent_proband_0.01CH_ConservationScore_Counts) <- SSC_parent_proband_0.01CH_ConservationScore_Counts$Group
SSC_parent_proband_0.01CH_ConservationScore_Counts <- SSC_parent_proband_0.01CH_ConservationScore_Counts[,-1]

#Fisher's exact test
SSC_parent_proband_0.01CH_ConservationScore_Fisher <- fisher.test(SSC_parent_proband_0.01CH_ConservationScore_Counts,
                                                        alternative='greater')

res <- rbind(res, data.frame(Name='SSC_parent_proband_0.01CH',
                             Child_Passed=SSC_parent_proband_0.01CH_ConservationScore_Counts$Passed[1],
                             Child_Failed=SSC_parent_proband_0.01CH_ConservationScore_Counts$Failed[1],
                             Parent_Passed=SSC_parent_proband_0.01CH_ConservationScore_Counts$Passed[2],
                             Parent_Failed=SSC_parent_proband_0.01CH_ConservationScore_Counts$Failed[2],
                             p_val=SSC_parent_proband_0.01CH_ConservationScore_Fisher$p.value,
                             OR=SSC_parent_proband_0.01CH_ConservationScore_Fisher$estimate))


SSC_parent_proband_0.1CH_ConservationScore_Groups <- Get_ConservationScore_Groups(SSC_meta,
                                                                                   SSC_parent_proband_0.1CH_SNVs,
                                                                                   SSC_proband_IDs,
                                                                                   SSC_father_IDs,
                                                                                   SSC_mother_IDs)

SSC_parent_proband_0.1CH_ConservationScore_Counts <- data.frame(SSC_parent_proband_0.1CH_ConservationScore_Groups %>% group_by(Group) %>%
                                                                   summarise(Passed=sum(Pass), Failed=sum(!Pass)))
rownames(SSC_parent_proband_0.1CH_ConservationScore_Counts) <- SSC_parent_proband_0.1CH_ConservationScore_Counts$Group
SSC_parent_proband_0.1CH_ConservationScore_Counts <- SSC_parent_proband_0.1CH_ConservationScore_Counts[,-1]

#Fisher's exact test
SSC_parent_proband_0.1CH_ConservationScore_Fisher <- fisher.test(SSC_parent_proband_0.1CH_ConservationScore_Counts,
                                                                  alternative='greater')

res <- rbind(res, data.frame(Name='SSC_parent_proband_0.1CH',
                             Child_Passed=SSC_parent_proband_0.1CH_ConservationScore_Counts$Passed[1],
                             Child_Failed=SSC_parent_proband_0.1CH_ConservationScore_Counts$Failed[1],
                             Parent_Passed=SSC_parent_proband_0.1CH_ConservationScore_Counts$Passed[2],
                             Parent_Failed=SSC_parent_proband_0.1CH_ConservationScore_Counts$Failed[2],
                             p_val=SSC_parent_proband_0.1CH_ConservationScore_Fisher$p.value,
                             OR=SSC_parent_proband_0.1CH_ConservationScore_Fisher$estimate))



SSC_parent_proband_CH_ConservationScore_Groups <- Get_ConservationScore_Groups(SSC_meta,
                                                                                  SSC_parent_proband_proc_SNVs,
                                                                                  SSC_proband_IDs,
                                                                                  SSC_father_IDs,
                                                                                  SSC_mother_IDs)

SSC_parent_proband_CH_ConservationScore_Counts <- data.frame(SSC_parent_proband_CH_ConservationScore_Groups %>% group_by(Group) %>%
                                                                  summarise(Passed=sum(Pass), Failed=sum(!Pass)))
rownames(SSC_parent_proband_CH_ConservationScore_Counts) <- SSC_parent_proband_CH_ConservationScore_Counts$Group
SSC_parent_proband_CH_ConservationScore_Counts <- SSC_parent_proband_CH_ConservationScore_Counts[,-1]

#Fisher's exact test
SSC_parent_proband_CH_ConservationScore_Fisher <- fisher.test(SSC_parent_proband_CH_ConservationScore_Counts,
                                                                 alternative='greater')

res <- rbind(res, data.frame(Name='SSC_parent_proband_CH',
                             Child_Passed=SSC_parent_proband_CH_ConservationScore_Counts$Passed[1],
                             Child_Failed=SSC_parent_proband_CH_ConservationScore_Counts$Failed[1],
                             Parent_Passed=SSC_parent_proband_CH_ConservationScore_Counts$Passed[2],
                             Parent_Failed=SSC_parent_proband_CH_ConservationScore_Counts$Failed[2],
                             p_val=SSC_parent_proband_CH_ConservationScore_Fisher$p.value,
                             OR=SSC_parent_proband_CH_ConservationScore_Fisher$estimate))


rownames(res) <- NULL
write.csv(res, '../DT/Conservation Score Plots/conservation_score_1cutoff_fishers_test.csv')
