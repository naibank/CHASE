library(yaml)
library(survival)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")


Combine_Tables <- function(table_names) {
  
  res <- data.frame()
  for (i in 1:length(table_names)) {
    res <- rbind(res, read.delim(table_names[[i]], stringsAsFactors = F))
  }
  
  return (res)
}

Get_Filtered_Metadata <- function(metadata_file, full_CNV_SNV_data_files, 
                                  child_relation='proband', comparison='parent-child') {
  child_list <- c()
  for (j in 1:length(full_CNV_SNV_data_files)) {
    message("--")
    message(j)
    data_dt <- read_yaml(full_CNV_SNV_data_files[[j]])
    for(i in 1:length(data_dt)){
      message(i)
      cnv <- data.frame(data_dt[[i]]$CNVs)
      if(!("Paternal" %in% cnv$Inheritance & "Maternal" %in% cnv$Inheritance)){
        child_list <- c(child_list, cnv$Sample.ID[1])
      }
    }
  }
  
  meta <- read.delim(metadata_file, stringsAsFactors = F)
  fam <- meta$Family.ID[meta$Sample.ID %in% child_list]
  if (comparison == 'parent-child') {
    meta <- meta[which(( meta$Sample.ID %in% child_list) |
                         (meta$Relation %in% c("father", "mother") & meta$Family.ID %in% fam)), ]
  }
  else if (comparison == 'child-child') {
    meta <- meta[which(meta$Sample.ID %in% child_list), ]
  }
  
  meta$Status <- ifelse(meta$Relation == child_relation, 1, 0)
  return (meta)
}

Get_Combined_SNVs <- function(parents_set, child_set) {
  parents <- Combine_Tables(parents_set)
  child <- Combine_Tables(child_set)
  
  snvs <- rbind(parents, child)
  
  return (snvs)
}

Fill_Remaining_CNV_Exonic_Sizes <- function(dt, cnvs) {
  for (i in 1:nrow(dt))
    if (is.na(dt$exonicSize[i]) & dt$'Sample.ID'[i] %in% cnvs$'Sample.ID') {
      dt$exonicSize[i] <- cnvs$exonicSize[which(cnvs$Sample.ID == dt$'Sample.ID'[i])]
    }
  
  return (dt)
}
Get_Logistic_Regression <- function(meta, snvs, gene_set_data, include_exon_size=F) {
  #Load additional CNV exonic size table
  if (include_exon_size == T) {
    CNVs_exonic_size <- rbind((read.delim("../DT/SSC_ProbandUNCNVs_exonicsizes.tsv",stringsAsFactors = F)),
                               (read.delim("../DT/SSC_UnaffectedSiblingsCNVs_exonicsizes.tsv",stringsAsFactors = F)))
    CNVs_exonic_size <- CNVs_exonic_size[!duplicated(CNVs_exonic_size[, 'Sample.ID']), ]
  }
  
  ### All variants
  all_var_count <- data.frame(table(snvs$X.Sample[!duplicated(snvs[, c("X.Sample", "X.id")])]))
  if (include_exon_size == T) {
    all_var_count$exonicSize <- snvs$exonicSize[match(all_var_count$Var1, snvs$X.Sample)]
  }
  
  dt.test <- merge(meta, all_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
  dt.test$Freq[is.na(dt.test$Freq)] <- 0
  
  if (include_exon_size == T) {
    dt.test <- Fill_Remaining_CNV_Exonic_Sizes(dt.test, CNVs_exonic_size)
    dt.test$exonicSize[is.na(dt.test$exonicSize)] <- 0
    all_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize, dt.test)
    all_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize + Freq, dt.test)
  }
  else {
    all_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
    all_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
  }

  all_ano <- anova(all_ref.lm, all_add.lm, method = "Chisq")
  
  ### Synonymous variants
  syn <- snvs[which(snvs$effect_priority == "synonymous SNV"), ]
  syn_var_count <- data.frame(table(syn$X.Sample[!duplicated(syn[, c("X.Sample", "X.id")])]))
  if (include_exon_size == T) {
    syn_var_count$exonicSize <- syn$exonicSize[match(syn_var_count$Var1, syn$X.Sample)]
  }
  
  dt.test <- merge(meta, syn_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
  dt.test$Freq[is.na(dt.test$Freq)] <- 0
  
  if (include_exon_size == T) {
    dt.test <- Fill_Remaining_CNV_Exonic_Sizes(dt.test, CNVs_exonic_size)
    dt.test$exonicSize[is.na(dt.test$exonicSize)] <- 0
    syn_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize, dt.test)
    syn_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize + Freq, dt.test)
  }
  else {
    syn_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
    syn_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
  }
  
  syn_ano <- anova(syn_ref.lm, syn_add.lm, method = "Chisq")
  
  ### Missense variants
  mis <- snvs[which(snvs$effect_priority == "nonsynonymous SNV"), ]
  mis_var_count <- data.frame(table(mis$X.Sample[!duplicated(mis[, c("X.Sample", "X.id")])]))
  if (include_exon_size == T) {
    mis_var_count$exonicSize <- mis$exonicSize[match(mis_var_count$Var1, mis$X.Sample)]
  }
  
  dt.test <- merge(meta, mis_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
  dt.test$Freq[is.na(dt.test$Freq)] <- 0
  
  if (include_exon_size == T) {
    dt.test <- Fill_Remaining_CNV_Exonic_Sizes(dt.test, CNVs_exonic_size)
    dt.test$exonicSize[is.na(dt.test$exonicSize)] <- 0
    mis_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize, dt.test)
    mis_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize + Freq, dt.test)
  }
  else {
    mis_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
    mis_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
  }
  
  mis_ano <- anova(mis_ref.lm, mis_add.lm, method = "Chisq")
  
  ### LOF variants
  lofs <- snvs[which(snvs$LoF), ]
  lof_var_count <- data.frame(table(lofs$X.Sample[!duplicated(lofs[, c("X.Sample", "X.id")])]))
  if (include_exon_size == T) {
    lof_var_count$exonicSize <- lofs$exonicSize[match(lof_var_count$Var1, lofs$X.Sample)]
  }
  
  dt.test <- merge(meta, lof_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
  dt.test$Freq[is.na(dt.test$Freq)] <- 0
  
  if (include_exon_size == T) {
    dt.test <- Fill_Remaining_CNV_Exonic_Sizes(dt.test, CNVs_exonic_size)
    dt.test$exonicSize[is.na(dt.test$exonicSize)] <- 0
    lof_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize, dt.test)
    lof_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize + Freq, dt.test)
  }
  else {
    lof_ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
    lof_add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
  }
  
  lof_ano <- anova(lof_ref.lm, lof_add.lm, method = "Chisq")
  
  main_variants_lm <- list(all=list(all_ref.lm, all_add.lm, all_ano),
                           syn=list(syn_ref.lm, syn_add.lm, syn_ano),
                           mis=list(mis_ref.lm, mis_add.lm, mis_ano),
                           lof=list(lof_ref.lm, lof_add.lm, lof_ano))
  
  ### Gene sets
  load(gene_set_data)
  gs_lm <- list()
  for (i in 1:length(gsMain)) {
    message(i)
    gs <- gsMain[[i]]
    gs_name <- names(gsMain)[i]
    gs_snvs <- snvs[which(snvs$entrez_id %in% gs), ]
    
    if (nrow(gs_snvs) > 0) {
      gs_var_count <- data.frame(table(gs_snvs$X.Sample[!duplicated(gs_snvs[, c("X.Sample", "X.id")])]))
      if (include_exon_size == T) {
        gs_var_count$exonicSize <- gs_snvs$exonicSize[match(gs_var_count$Var1, gs_snvs$X.Sample)]
      }
      
      dt.test <- merge(meta, gs_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
      dt.test$Freq[is.na(dt.test$Freq)] <- 0
      
      if (include_exon_size == T) {
        dt.test <- Fill_Remaining_CNV_Exonic_Sizes(dt.test, CNVs_exonic_size)
        dt.test$exonicSize[is.na(dt.test$exonicSize)] <- 0
        ref.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize, dt.test)
        add.lm <- clogit(Status ~ strata(Family.ID) + Sex + exonicSize + Freq, dt.test)
      }
      else {
        ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
        add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
      }
      
      ano <- anova(ref.lm, add.lm, method = "Chisq")
      
      gs_lm[[gs_name]] <- list(ref.lm, add.lm, ano)
    }
  }
  
  return (list(main_variants_lm, gs_lm))
}

gene_set_data_path = "gsMain_PGC_2021.RData"
### MSSNG Parent-Proband
MSSNG_metadata_path <- "MSSNG_metadata.tsv"
MSSNG_full_data_files <- list("MSSNG_ILMN_CH_Data_CNV10P_SNV.yaml", "MSSNG_CG_CH_Data_CNV10P_SNV.yaml")
MSSNG_meta <- Get_Filtered_Metadata(MSSNG_metadata_path, MSSNG_full_data_files, 
                                    child_relation = 'proband')

MSSNG_parent_file_list <- list("../DT/ILMN_FatherSNVs_in_PaternalCNVs.tsv", "../DT/CG_FatherSNVs_in_PaternalCNVs.tsv",
                         "../DT/ILMN_MotherSNVs_in_MaternalCNVs.tsv", "../DT/CG_MotherSNVs_in_MaternalCNVs.tsv")
MSSNG_proband_file_list <- list("../DT/ILMN_ProbandSNVs_in_PaternalCNVs.tsv", "../DT/CG_ProbandSNVs_in_PaternalCNVs.tsv",
                                "../DT/ILMN_ProbandSNVs_in_MaternalCNVs.tsv", "../DT/CG_ProbandSNVs_in_MaternalCNVs.tsv")
MSSNG_parent_proband_snvs <- Get_Combined_SNVs(MSSNG_parent_file_list, MSSNG_proband_file_list)
MSSNG_parent_proband_regression_res <- Get_Logistic_Regression(MSSNG_meta, 
                                                               MSSNG_parent_proband_snvs,
                                                               gene_set_data = gene_set_data_path)


### SSC Parent-Proband
SSC_metadata_path <- "SSC_metadata.tsv"
SSC_full_data_files <- list("SSC_CH_Data_CNV10P_SNV.yaml")
SSC_meta <- Get_Filtered_Metadata(SSC_metadata_path, SSC_full_data_files, 
                                    child_relation = 'proband')

SSC_parent_of_proband_file_list <- list("../DT/SSC_FatherSNVs_in_PaternalCNVs.tsv", "../DT/SSC_MotherSNVs_in_MaternalCNVs.tsv")
SSC_proband_file_list <- list("../DT/SSC_ProbandSNVs_in_PaternalCNVs.tsv", "../DT/SSC_ProbandSNVs_in_MaternalCNVs.tsv")
SSC_parent_proband_snvs <- Get_Combined_SNVs(SSC_parent_of_proband_file_list, 
                                             SSC_proband_file_list)
SSC_parent_proband_regression_res <- Get_Logistic_Regression(SSC_meta, 
                                                             SSC_parent_proband_snvs,
                                                             gene_set_data = gene_set_data_path)

### SSC Parent-Unaffected Siblings
SSC_US_metadata_path <- "SSC_metadata.tsv"
SSC_US_full_data_files <- list("SSC_CH.unaffectedSiblings_Data_CNV10P_SNV.yaml")
SSC_US_meta <- Get_Filtered_Metadata(SSC_US_metadata_path, SSC_US_full_data_files, 
                                  child_relation = 'unaffected sibling')

SSC_parent_of_unaffectedSibling_file_list <- list("../DT/SSC_FatherUSSNVs_in_PaternalCNVs.tsv", "../DT/SSC_MotherUSSNVs_in_MaternalCNVs.tsv")
SSC_US_file_list <- list("../DT/SSC_unaffectedSibSNVs_in_PaternalCNVs.tsv", "../DT/SSC_unaffectedSibSNVs_in_MaternalCNVs.tsv")
SSC_parent_US_snvs <- Get_Combined_SNVs(SSC_parent_of_unaffectedSibling_file_list, 
                                        SSC_US_file_list)
SSC_parent_US_regression_res <- Get_Logistic_Regression(SSC_US_meta, 
                                                        SSC_parent_US_snvs,
                                                        gene_set_data = gene_set_data_path)

### SSC Proband-Unaffected Siblings
SSC_sibilng_metadata_path <- "SSC_metadata.tsv"
SSC_sibling_full_data_files <- list("SSC_CH_Data_CNV10P_SNV.yaml","SSC_CH.unaffectedSiblings_Data_CNV10P_SNV.yaml")
SSC_sibling_meta <- Get_Filtered_Metadata(SSC_sibilng_metadata_path, SSC_sibling_full_data_files, 
                                     child_relation = 'proband',
                                     comparison = 'child-child')

SSC_probandSibling_file_list <- list("../DT/SSC_ProbandUNSNVs.tsv")
SSC_USSibilng_file_list <- list("../DT/SSC_UnaffectedSiblingsSNVs.tsv")
SSC_proband_US_snvs <- Get_Combined_SNVs(SSC_probandSibling_file_list, 
                                         SSC_USSibilng_file_list)
SSC_proband_US_regression_res <- Get_Logistic_Regression(SSC_sibling_meta, 
                                                         SSC_proband_US_snvs,
                                                         gene_set_data = gene_set_data_path,
                                                         include_exon_size = T)


#### CH Event Limited to 0.1% *NOTE: Overwrites previous vars*####
### MSSNG Parent-Proband
MSSNG_parent_file_list <- list("../DT/ILMN_FatherSNVs_in_PaternalCNVs_01CHfreq.tsv", "../DT/CG_FatherSNVs_in_PaternalCNVs_01CHfreq.tsv",
                               "../DT/ILMN_MotherSNVs_in_MaternalCNVs_01CHfreq.tsv", "../DT/CG_MotherSNVs_in_MaternalCNVs_01CHfreq.tsv")
MSSNG_proband_file_list <- list("../DT/ILMN_ProbandSNVs_in_PaternalCNVs_01CHfreq.tsv", "../DT/CG_ProbandSNVs_in_PaternalCNVs_01CHfreq.tsv",
                                "../DT/ILMN_ProbandSNVs_in_MaternalCNVs_01CHfreq.tsv", "../DT/CG_ProbandSNVs_in_MaternalCNVs_01CHfreq.tsv")
MSSNG_parent_proband_snvs <- Get_Combined_SNVs(MSSNG_parent_file_list, MSSNG_proband_file_list)
MSSNG_parent_proband_regression_res <- Get_Logistic_Regression(MSSNG_meta, 
                                                               MSSNG_parent_proband_snvs,
                                                               gene_set_data = gene_set_data_path)


### SSC Parent-Proband
SSC_parent_of_proband_file_list <- list("../DT/SSC_FatherSNVs_in_PaternalCNVs_01CHfreq.tsv", "../DT/SSC_MotherSNVs_in_MaternalCNVs_01CHfreq.tsv")
SSC_proband_file_list <- list("../DT/SSC_ProbandSNVs_in_PaternalCNVs_01CHfreq.tsv", "../DT/SSC_ProbandSNVs_in_MaternalCNVs_01CHfreq.tsv")
SSC_parent_proband_snvs <- Get_Combined_SNVs(SSC_parent_of_proband_file_list, 
                                             SSC_proband_file_list)
SSC_parent_proband_regression_res <- Get_Logistic_Regression(SSC_meta, 
                                                             SSC_parent_proband_snvs,
                                                             gene_set_data = gene_set_data_path)

### SSC Parent-Unaffected Siblings
SSC_parent_of_unaffectedSibling_file_list <- list("../DT/SSC_FatherUSSNVs_in_PaternalCNVs_01CHfreq.tsv", "../DT/SSC_MotherUSSNVs_in_MaternalCNVs_01CHfreq.tsv")
SSC_US_file_list <- list("../DT/SSC_unaffectedSibSNVs_in_PaternalCNVs_01CHfreq.tsv", "../DT/SSC_unaffectedSibSNVs_in_MaternalCNVs_01CHfreq.tsv")
SSC_parent_US_snvs <- Get_Combined_SNVs(SSC_parent_of_unaffectedSibling_file_list, 
                                        SSC_US_file_list)
SSC_parent_US_regression_res <- Get_Logistic_Regression(SSC_US_meta, 
                                                        SSC_parent_US_snvs,
                                                        gene_set_data = gene_set_data_path)

### SSC Proband-Unaffected Siblings
SSC_probandSibling_file_list <- list("../DT/SSC_ProbandUNSNVs_01CHfreq.tsv")
SSC_USSibilng_file_list <- list("../DT/SSC_UnaffectedSiblingsSNVs_01CHfreq.tsv")
SSC_proband_US_snvs <- Get_Combined_SNVs(SSC_probandSibling_file_list, 
                                         SSC_USSibilng_file_list)
SSC_proband_US_regression_res <- Get_Logistic_Regression(SSC_sibling_meta, 
                                                         SSC_proband_US_snvs,
                                                         gene_set_data = gene_set_data_path,
                                                         include_exon_size = T)



Generate_Gene_Set_Dataframe <- function(regression, coeff=2) {
  gene_sets_list <- regression[[2]]
  gene_sets_names <- names(gene_sets_list)
  res <- data.frame()
  for (i in 1:length(gene_sets_list)) {
    name <- gene_sets_names[[i]]
    gs <- gene_sets_list[[i]]
    row <- data.frame(GeneSet=name, Coeff=gs[[2]]$coefficients[coeff], PVal=gs[[3]]$`P(>|Chi|)`[2])
    res <- rbind(res, row)
  }
  
  res <- res[order(res$PVal),]
  res <- res[1:10, ]
  
  return (res)
}

MSSNG_parent_proband_gs_df <- Generate_Gene_Set_Dataframe(MSSNG_parent_proband_regression_res)
SSC_parent_proband_gs_df <- Generate_Gene_Set_Dataframe(SSC_parent_proband_regression_res)
SSC_parent_US_gs_df <- Generate_Gene_Set_Dataframe(SSC_parent_US_regression_res)
SSC_proband_US_gs_df <- Generate_Gene_Set_Dataframe(SSC_proband_US_regression_res, coeff=3)

# proband <- rbind(read.delim("ILMN_ProbandSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F),
#                  read.delim("ILMN_ProbandSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F))
# 
# parent <- rbind(read.delim("ILMN_FatherSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F),
#                 read.delim("ILMN_MotherSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F))
# 
# yaml_dt <- read_yaml("MSSNG_ILMN_CH_Data_CNV10P_SNV.yaml")
# 
# proband_list <- c()
# for(i in 1:length(yaml_dt)){
#   cnv <- data.frame(yaml_dt[[i]]$CNVs)
#   if(!("Paternal" %in% cnv$Inheritance & "Maternal" %in% cnv$Inheritance)){
#     proband_list <- c(proband_list, cnv$Sample.ID[1])
#   }
# }
# 
# meta <- read.delim("MSSNG_metadata.tsv", stringsAsFactors = F)
# fam <- meta$Family.ID[meta$Sample.ID %in% proband_list]
# 
# 
# meta <- meta[which(( meta$Sample.ID %in% proband_list) |
#                (meta$Relation %in% c("father", "mother") & meta$Family.ID %in% fam)), ]
# meta$Status <- ifelse(meta$Relation == "proband", 1, 0)
# 
# snvs <- rbind(proband, parent)
# 
# ### all variants
# all_var_count <- data.frame(table(snvs$X.Sample[!duplicated(snvs[, c("X.Sample", "X.id")])]))
# 
# dt.test <- merge(meta, all_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
# dt.test$Freq[is.na(dt.test$Freq)] <- 0
# 
# ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
# add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
# 
# ano <- anova(ref.lm, add.lm, method = "Chisq")
# 
# ###coefficient
# add.lm$coefficients[2]
# 
# ###OR
# exp(add.lm$coefficients[2])
# 
# ### p-value
# ano$`P(>|Chi|)`[2]
# 
# 
# ### lof variants
# lofs <- snvs[which(snvs$LoF), ]
# lof_var_count <- data.frame(table(lofs$X.Sample[!duplicated(lofs[, c("X.Sample", "X.id")])]))
# 
# dt.test <- merge(meta, lof_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
# dt.test$Freq[is.na(dt.test$Freq)] <- 0
# 
# ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
# add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
# 
# ano <- anova(ref.lm, add.lm, method = "Chisq")
# 
# ###coeffcient
# add.lm$coefficients[2]
# 
# ### p-value
# ano$`P(>|Chi|)`[2]
# 
# ##### gene set analysis
# load("gsMain_PGC_2021.RData")
# names(gsMain)
# 
# gs <- gsMain$Neurof_UnionInclusive
# gs_snvs <- snvs[which(snvs$entrez_id %in% gs), ]
# 
# gs_var_count <- data.frame(table(gs_snvs$X.Sample[!duplicated(gs_snvs[, c("X.Sample", "X.id")])]))
# 
# dt.test <- merge(meta, gs_var_count, by.x = "Sample.ID", by.y = "Var1", all.x = T)
# dt.test$Freq[is.na(dt.test$Freq)] <- 0
# 
# ref.lm <- clogit(Status ~ strata(Family.ID) + Sex, dt.test)
# add.lm <- clogit(Status ~ strata(Family.ID) + Sex + Freq, dt.test)
# 
# ano <- anova(ref.lm, add.lm, method = "Chisq")
# 
# ###coeffcient
# add.lm$coefficients[2]
# 
# ### p-value
# ano$`P(>|Chi|)`[2]
