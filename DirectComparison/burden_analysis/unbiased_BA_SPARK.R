library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

Determine_CNV_MAF <- function(CNVs) {
  CNV_MAF <- pmax(CNVs$CGparentalPercFreq_50percRecipOverlap,
                  pmax(CNVs$cnvnIlmXParentalPercFreq_50percRecipOverlap,
                       pmax(CNVs$erdsIlmXParentalPercFreq_50percRecipOverlap,
                            pmax(CNVs$cnvnIlm2ParentalPercFreq_50percRecipOverlap,
                                 pmax(CNVs$erdsIlm2ParentalPercFreq_50percRecipOverlap,
                                      pmax(CNVs$otgCnvnPercFreq_50percRecipOverlap,
                                           pmax(CNVs$otgErdsPercFreq_50percRecipOverlap,
                                                pmax(CNVs$sscErdsPercFreq_50percRecipOverlap,
                                                     CNVs$sscCnvnPercFreq_50percRecipOverlap, na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T)/100
  return (CNV_MAF)
}


Get_Failed_QC_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("MSSNG+SSC.CNVs.tsv")
  ls_exclude <- unique(rbind(ls1, ls2))
  
  full_metadata <- fread(file, data.table = F)
  familyID_exclude <- full_metadata$`Family ID`[which(full_metadata$`Sample ID` %in% ls_exclude$Sample)]
  sampleID_exclude <- full_metadata$`Sample ID`[which(full_metadata$`Family ID` %in% familyID_exclude)]
  
  return (sampleID_exclude)
}

Get_Filtered_Metadata <- function(metadata_file, child_relation='proband', 
                                  comparison='parent-child') {
  meta <- read.delim(metadata_file, stringsAsFactors = F)
  
  # Get samples that failed QC
  failed_QC <- Get_Failed_QC_Samples()
  # Get samples with cr CNV
  samples_with_cr_CNV <- Get_cr_CNV_Samples()
  # Get families to exclude that have parents AS the proband as well
  families_to_exclude <- meta$Family.ID[(meta$Sample.ID[meta$Relation == 'proband'] %in% meta$Mother.ID) |
                                          (meta$Sample.ID[meta$Relation == 'proband'] %in% meta$Father.ID)]
  # Get samples with ASD associated CNVs and LOFs
  ASD_associated_samples <- Get_SampleID_ASD_Associated_CNVs_and_LOFs(metadata_file)
  
  # Filter out from metadata
  meta <- meta[!meta$Sample.ID %in% failed_QC, ]
  meta <- meta[!meta$Sample.ID %in% samples_with_cr_CNV, ]
  meta <- meta[!meta$Family.ID %in% families_to_exclude, ]
  meta <- meta[!meta$Sample.ID %in% ASD_associated_samples, ]
  
  child_relations <- ifelse(child_relation == 'proband', c('proband'), c('unaffected sibling', 'other sibling'))
  # Filter out irrelevant child if parent-child comparison
  if (comparison == 'parent-child') {
    meta <- meta[meta$Relation %in% child_relations | meta$Relation == 'mother'
                 | meta$Relation == 'father', ]
    # Get child-parent trios only
    IDs_relations <- meta[, c('Family.ID', 'Sample.ID', 'Relation')]
    IDs_relations$FamRel <- paste(IDs_relations$Family.ID, IDs_relations$Relation, sep='.')
    unique_famIDs_and_relations <- dplyr::distinct(IDs_relations,FamRel,
                                                   .keep_all=T)
    famID_freqs <- data.frame(table(unique_famIDs_and_relations$Family.ID))
    trios_IDs <- unique_famIDs_and_relations$Sample.ID[unique_famIDs_and_relations$Family.ID %in%
                                                         famID_freqs$Var1[famID_freqs$Freq == 3]] # 1 child, 1 father, 1 mother
    meta <- meta[meta$Sample.ID %in% trios_IDs, ]
  }
  # Filter out parents if child-child comparison
  else if (comparison == 'child-child') {
    meta <- meta[meta$Relation %in% c('proband', 'unaffected sibling', 'other sibling'), ]
    # Get only samples from families with 1 of each sibling
    proband_fam_IDs <- meta$Family.ID[meta$Relation == 'proband']
    sibling_fam_IDs <- meta$Family.ID[meta$Relation %in% c('unaffected sibling', 'other sibling')]
    meta <- meta[meta$Family.ID %in% proband_fam_IDs &
                   meta$Family.ID %in% sibling_fam_IDs,]
  }
  
  meta$Status <- ifelse(meta$Relation %in% child_relations, 1, 0)
  return (meta)
}

Get_Clogit_ANOVA <- function(meta, snvs, ref, mod, type='clogit') {
  dt.test <- merge(meta, snvs, by.x = "Sample.ID", by.y = "Var1", all.x = T)
  dt.test$Freq[is.na(dt.test$Freq)] <- 0
  dt.test$exonicSize[is.na(dt.test$exonicSize)] <- 0
  dt.test$totalCHEvents[is.na(dt.test$totalCHEvents)] <- 0 # TODO: Change this to reset every NA cell instead
  dt.test$totalCHEvents_syn[is.na(dt.test$totalCHEvents_syn)] <- 0 # TODO: Change this to reset every NA cell instead
  dt.test$CHrate <- 0
  dt.test$CHrate[!dt.test$Freq == 0] <- dt.test$Freq[!dt.test$Freq == 0] / dt.test$exonicSize[!dt.test$Freq == 0]
  
  if (type == 'clogit') {
    ref.lm <- clogit(ref, dt.test)
    add.lm <- clogit(mod, dt.test)
  }
  else if (type == 'glm') {
    ref.lm <- glm(ref, dt.test, family='binomial')
    add.lm <- glm(mod, dt.test, family='binomial')
  }
  
  ano <- anova(ref.lm, add.lm, test = "Chisq")
  
  return (list(ref.lm, add.lm, ano))
}

Get_Log_Reg <- function(meta, snvs, gene_set_data, exonic_size_data,
                        ref = Status ~ Sex, mod = Status ~ Sex + CHrate,
                        snv_count_groups = 'parent-child', type='clogit',
                        apply_SNV_restrictions=T) {
  snvs <- snvs[snvs$X.Sample %in% meta$Sample.ID,]
  if (snv_count_groups == 'parent-child') {
    group1_IDs <- meta$Sample.ID[!meta$Relation == 'mother' & !meta$Relation == 'father']
    group2_IDs <- meta$Sample.ID[meta$Relation == 'father']
    group3_IDs <- meta$Sample.ID[meta$Relation == 'mother']
  }
  else if (snv_count_groups == 'proband-unaffected sibling') {
    group1_IDs <- meta$Sample.ID[meta$Relation == 'proband']
    group2_IDs <- meta$Sample.ID[meta$Relation == 'unaffected sibling']
    group3_IDs <- meta$Sample.ID[meta$Relation == 'unaffected sibling'] # placeholder
  }
  
  all_var <- data.frame(table(snvs$X.Sample[!duplicated(snvs[, c("X.Sample", "X.id")])]))
  all_syn_var <- all_var[all_var$effect_priority == "synonymous SNV", c('Var1', 'Freq')]
  all_var <- all_var[, c('Var1', 'Freq')]
  colnames(all_var)[colnames(all_var)=='Freq'] <- 'totalCHEvents'
  colnames(all_syn_var)[colnames(all_syn_var)=='Freq'] <- 'totalCHEvents_syn'
  
  all_reg_base <- list() # Placeholder for base all variant reg with no restrictions
  if (apply_SNV_restrictions) { # Apply to each type alone by separating all variant counts first
    SNV_count <- data.frame(Group1 = length(snvs$X.Sample[snvs$X.Sample %in% group1_IDs]),
                            Group2 = length(snvs$X.Sample[snvs$X.Sample %in% group2_IDs]),
                            Group3 = length(snvs$X.Sample[snvs$X.Sample %in% group3_IDs]))
    all_var_count_base <- data.frame(table(snvs$X.Sample[!duplicated(snvs[, c("X.Sample", "X.id")])]))
    all_var_count_base <- merge(all_var_count_base, exonic_size_data, by.x='Var1', by.y='Sample.ID', all.x=T)
    all_var_count_base <- merge(all_var_count_base, all_var, by.x='Var1', by.y='Var1', all=T)
    all_var_count_base <- merge(all_var_count_base, all_syn_var, by.x='Var1', by.y='Var1', all=T)
    all_var_count_base[is.na(all_var_count_base)] <- 0
    all_reg_base <- Get_Clogit_ANOVA(meta, all_var_count_base, ref, mod, type)
    all_reg_base <- append(all_reg_base, SNV_count)
    message('Done base all var reg')
    
    snvs <- snvs[snvs$gnomAD_pRec >= 0.9 & snvs$gnomAD_oe_lof_upper >= 0.35, ]
  }
  
  ### All variants
  SNV_count <- data.frame(Group1 = length(snvs$X.Sample[snvs$X.Sample %in% group1_IDs]),
                          Group2 = length(snvs$X.Sample[snvs$X.Sample %in% group2_IDs]),
                          Group3 = length(snvs$X.Sample[snvs$X.Sample %in% group3_IDs]))
  all_var_count <- data.frame(table(snvs$X.Sample[!duplicated(snvs[, c("X.Sample", "X.id")])]))
  all_var_count <- merge(all_var_count, exonic_size_data, by.x='Var1', by.y='Sample.ID', all.x=T)
  all_var_count <- merge(all_var_count, all_var, by.x='Var1', by.y='Var1', all=T)
  all_var_count <- merge(all_var_count, all_syn_var, by.x='Var1', by.y='Var1', all=T)
  all_var_count[is.na(all_var_count)] <- 0
  all_reg <- Get_Clogit_ANOVA(meta, all_var_count, ref, mod, type)
  all_reg <- append(all_reg, SNV_count)
  message("Done all variant")
  ### Synonymous variants
  syn <- snvs[which(snvs$effect_priority == "synonymous SNV"), ]
  syn_reg <- data.frame()
  if (nrow(syn) > 0) {
    SNV_count <- data.frame(Group1 = length(syn$X.Sample[syn$X.Sample %in% group1_IDs]),
                            Group2 = length(syn$X.Sample[syn$X.Sample %in% group2_IDs]),
                            Group3 = length(syn$X.Sample[syn$X.Sample %in% group3_IDs]))
    syn_var_count <- data.frame(table(syn$X.Sample[!duplicated(syn[, c("X.Sample", "X.id")])]))
    syn_var_count <- merge(syn_var_count, exonic_size_data, by.x='Var1', by.y='Sample.ID', all.x=T)
    syn_var_count <- merge(syn_var_count, all_var, by.x='Var1', by.y='Var1', all=T)
    syn_var_count <- merge(syn_var_count, all_syn_var, by.x='Var1', by.y='Var1', all=T)
    syn_var_count[is.na(syn_var_count)] <- 0
    syn_reg <- Get_Clogit_ANOVA(meta, syn_var_count, ref, mod, type)
    syn_reg <- append(syn_reg, SNV_count)
  }
  message("Done syn variant")
  # ### Missense variants
  mis <- snvs[which(snvs$effect_priority == "nonsynonymous SNV"), ]
  SNV_count <- data.frame(Group1 = length(mis$X.Sample[mis$X.Sample %in% group1_IDs]),
                          Group2 = length(mis$X.Sample[mis$X.Sample %in% group2_IDs]),
                          Group3 = length(mis$X.Sample[mis$X.Sample %in% group3_IDs]))
  mis_var_count <- data.frame(table(mis$X.Sample[!duplicated(mis[, c("X.Sample", "X.id")])]))
  mis_var_count <- merge(mis_var_count, exonic_size_data, by.x='Var1', by.y='Sample.ID', all.x=T)
  mis_var_count <- merge(mis_var_count, all_var, by.x='Var1', by.y='Var1', all=T)
  mis_var_count <- merge(mis_var_count, all_syn_var, by.x='Var1', by.y='Var1', all=T)
  mis_var_count[is.na(mis_var_count)] <- 0
  mis_reg <- Get_Clogit_ANOVA(meta, mis_var_count, ref, mod, type)
  mis_reg <- append(mis_reg, SNV_count)
  message("Done mis variant")
  # ### LOF variants
  lofs <- snvs[which(snvs$LoF), ]
  SNV_count <- data.frame(Group1 = length(lofs$X.Sample[lofs$X.Sample %in% group1_IDs]),
                          Group2 = length(lofs$X.Sample[lofs$X.Sample %in% group2_IDs]),
                          Group3 = length(lofs$X.Sample[lofs$X.Sample %in% group3_IDs]))
  lof_var_count <- data.frame(table(lofs$X.Sample[!duplicated(lofs[, c("X.Sample", "X.id")])]))
  lof_reg <- data.frame()
  if (nrow(lof_var_count) > 0) {
    lof_var_count <- merge(lof_var_count, exonic_size_data, by.x='Var1', by.y='Sample.ID', all.x=T)
    lof_var_count <- merge(lof_var_count, all_var, by.x='Var1', by.y='Var1', all=T)
    lof_var_count <- merge(lof_var_count, all_syn_var, by.x='Var1', by.y='Var1', all=T)
    lof_var_count[is.na(lof_var_count)] <- 0
    lof_reg <- Get_Clogit_ANOVA(meta, lof_var_count, ref, mod, type)
    lof_reg <- append(lof_reg, SNV_count)
  }
  message("Done lof variant")
  main_variants_lm <- list(all=all_reg,
                           syn=syn_reg,
                           mis=mis_reg,
                           lof=lof_reg,
                           all_base=all_reg_base)
  message("Done main variant regression")
  # ### Gene sets
  load(gene_set_data)
  gs_lm <- list()
  for (i in 1:length(gsMain)) {
    tryCatch(
      expr = {
        message(i)
        gs <- gsMain[[i]]
        gs_name <- names(gsMain)[i]
        gs_snvs <- snvs[which(snvs$entrez_id %in% gs &
                                (snvs$LoF | snvs$effect_priority == "nonsynonymous SNV")), ]
        if (nrow(gs_snvs) > 0) {
          SNV_count <- data.frame(Group1 = length(gs_snvs$X.Sample[gs_snvs$X.Sample %in% group1_IDs]),
                                  Group2 = length(gs_snvs$X.Sample[gs_snvs$X.Sample %in% group2_IDs]),
                                  Group3 = length(gs_snvs$X.Sample[gs_snvs$X.Sample %in% group3_IDs]))
          gs_var_count <- data.frame(table(gs_snvs$X.Sample[!duplicated(gs_snvs[, c("X.Sample", "X.id")])]))
          gs_var_count <- merge(gs_var_count, exonic_size_data, by.x='Var1', by.y='Sample.ID', all.x=T)
          gs_var_count <- merge(gs_var_count, all_var, by.x='Var1', by.y='Var1', all=T)
          gs_var_count <- merge(gs_var_count, all_syn_var, by.x='Var1', by.y='Var1', all=T)
          gs_var_count[is.na(gs_var_count)] <- 0
          gs_res <- Get_Clogit_ANOVA(meta, gs_var_count, ref, mod, type)
          gs_res <- append(gs_res, SNV_count)
          gs_lm[[gs_name]] <- gs_res
        }
      },
      error = function(e){
        message(paste('Skipping gene set: ', gs_name))
        message(e)
      },
      finally = {
        next
      }
    )
  }
  
  return (list(main_variants_lm, gs_lm))
}

gene_exon_data <- data.table::fread("hg38_refGene_20200708.exon.txt", data.table = F)
gene_set_data_path <- "gsMain_PGC_2021.RData"
### SPARK Parent-Proband 
SPARK_meta <- data.frame()
SPARK_parent_proband_SNVs <- data.frame()
SPARK_parent_proband_ExonicSizes <- data.frame()
SPARK_parent_proband_CNVs <- data.frame()

for (i in 1:3) {
  SPARK_metadata_path <- paste("./SPARK/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  SPARK_meta <- rbind(SPARK_meta, Get_Filtered_Metadata(SPARK_metadata_path, child_relation = 'proband'))
  
  SPARK_SNVsExonicSizes_path <- paste('./SPARK/SPARKWGS',i,'_parent_proband_SNVsExonicSizes.yaml',sep='')
  SPARK_parent_proband_SNVsExonicSizes_curr <- yaml::yaml.load_file(SPARK_SNVsExonicSizes_path)
  
  SPARK_parent_proband_SNVs <- rbind(SPARK_parent_proband_SNVs,
                                     data.frame(SPARK_parent_proband_SNVsExonicSizes_curr[[1]]))
  SPARK_parent_proband_ExonicSizes <- rbind(SPARK_parent_proband_ExonicSizes,
                                            data.frame(SPARK_parent_proband_SNVsExonicSizes_curr[[2]]))
  SPARK_parent_proband_CNVs <- rbind(SPARK_parent_proband_CNVs,
                                     data.frame(SPARK_parent_proband_SNVsExonicSizes_curr[[3]]))
}
SPARK_parent_proband_CNVs$freq_max <- Determine_CNV_MAF(SPARK_parent_proband_CNVs)

SPARK_parent_proband_SNVs <- SPARK_parent_proband_SNVs[SPARK_parent_proband_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                       | SPARK_parent_proband_SNVs$LoF,] 
SPARK_parent_proband_SNVs$UID <- paste(SPARK_parent_proband_SNVs$X.Sample, SPARK_parent_proband_SNVs$X.id, sep='.')
SPARK_parent_proband_SNVs <- SPARK_parent_proband_SNVs[!(SPARK_parent_proband_SNVs$LoF & SPARK_parent_proband_SNVs$freq_max > 0.01), ]
SPARK_parent_proband_SNVs <- SPARK_parent_proband_SNVs[!is.na(SPARK_parent_proband_SNVs$X.Sample),]

# MSSNG_parent_proband_SV_SNVsExonicSizes <- yaml::yaml.load_file("./SV/MSSNG_parent_proband_SV_SNVsExonicSizes.yaml")
# MSSNG_parent_proband_SV_SNVs <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[1]])
# MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[MSSNG_parent_proband_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
#                                                              | MSSNG_parent_proband_SV_SNVs$LoF,] 
# MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[MSSNG_parent_proband_SV_SNVs$cnvENDAnn - MSSNG_parent_proband_SV_SNVs$cnvSTARTAnn <= 10000, ]
# MSSNG_parent_proband_SV_SNVs$UID <- paste(MSSNG_parent_proband_SV_SNVs$X.Sample, MSSNG_parent_proband_SV_SNVs$X.id, sep='.')
# MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[!(MSSNG_parent_proband_SV_SNVs$LoF & MSSNG_parent_proband_SV_SNVs$freq_max > 0.01), ]
# MSSNG_parent_proband_SV_ExonicSizes <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[2]])

# #Combine
# MSSNG_parent_proband_all_ExonicSizes <- merge(MSSNG_parent_proband_ExonicSizes, MSSNG_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
# MSSNG_parent_proband_all_ExonicSizes$exonicSize.x[is.na(MSSNG_parent_proband_all_ExonicSizes$exonicSize.x)] <- 0
# MSSNG_parent_proband_all_ExonicSizes$exonicSize.y[is.na(MSSNG_parent_proband_all_ExonicSizes$exonicSize.y)] <- 0
# MSSNG_parent_proband_all_ExonicSizes$exonicSize <- pmax(MSSNG_parent_proband_all_ExonicSizes$exonicSize.x, MSSNG_parent_proband_all_ExonicSizes$exonicSize.y)
# MSSNG_parent_proband_combined_SNVs <- rbind(MSSNG_parent_proband_SNVs, MSSNG_parent_proband_SV_SNVs)
# MSSNG_parent_proband_proc_SNVs <- dplyr::distinct(MSSNG_parent_proband_combined_SNVs, UID, .keep_all=T)

# # Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
# MSSNG_parent_proband_SVs_ILMN <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[3]])
# MSSNG_parent_proband_SVs_ILMN <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$length <= 10000
#                                                                & MSSNG_parent_proband_SVs_ILMN$Sample.ID %in% MSSNG_meta$Sample.ID, ]
# MSSNG_parent_proband_SVs_CG <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[4]])
# MSSNG_parent_proband_SVs_CG <- MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$length <= 10000
#                                                            & MSSNG_parent_proband_SVs_CG$Sample.ID %in% MSSNG_meta$Sample.ID, ]
# MSSNG_homozyg_SVs <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
# MSSNG_homozyg_SVs <- rbind(MSSNG_homozyg_SVs,MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')])
# MSSNG_homozyg_SVs$SVUID <- with(MSSNG_homozyg_SVs, paste0(sample,CHROM,START,END))
# 
# MSSNG_parent_proband_proc_SNVs <- MSSNG_parent_proband_proc_SNVs[!with(MSSNG_parent_proband_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% MSSNG_homozyg_SVs$SVUID,]

SPARK_parent_proband_proc_SNVs <- SPARK_parent_proband_SNVs  # Temp until we get SPARK SVs and do real processing (like above)

SPARK_parent_proband_proc_SNVs$Min_Dist <- by(SPARK_parent_proband_proc_SNVs, seq_len(nrow(SPARK_parent_proband_proc_SNVs)),
                                              function(r) r$MIN = min(abs(r$POS - SPARK_parent_proband_proc_SNVs$POS[SPARK_parent_proband_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                                       SPARK_parent_proband_proc_SNVs$CHROM == r$CHROM &
                                                                                                                       SPARK_parent_proband_proc_SNVs$POS != r$POS])) > 50)

SPARK_parent_proband_proc_SNVs <- SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$Min_Dist,]

SPARK_parent_proband_clogit_res <- Get_Log_Reg(SPARK_meta,
                                               SPARK_parent_proband_proc_SNVs,
                                               gene_set_data_path,
                                               SPARK_parent_proband_ExonicSizes,
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

### SPARK Parent-Unaffected Siblings
SPARK_US_meta <- data.frame()
SPARK_parent_US_SNVs <- data.frame()
SPARK_parent_US_ExonicSizes <- data.frame()
SPARK_parent_US_CNVs <- data.frame()
for (i in 1:3) {
  SPARK_US_metadata_path <- paste("./SPARK/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  SPARK_US_meta <- rbind(SPARK_US_meta, Get_Filtered_Metadata(SPARK_US_metadata_path, child_relation = 'unaffected sibling'))
  
  SPARK_US_SNVsExonicSizes_path <- paste('./SPARK/SPARKWGS',i,'_parent_US_SNVsExonicSizes.yaml',sep='')
  SPARK_parent_US_SNVsExonicSizes_curr <- yaml::yaml.load_file(SPARK_US_SNVsExonicSizes_path)
  
  SPARK_parent_US_SNVs <- rbind(SPARK_parent_US_SNVs,
                                     data.frame(SPARK_parent_US_SNVsExonicSizes_curr[[1]]))
  SPARK_parent_US_ExonicSizes <- rbind(SPARK_parent_US_ExonicSizes,
                                            data.frame(SPARK_parent_US_SNVsExonicSizes_curr[[2]]))
  SPARK_parent_US_CNVs <- rbind(SPARK_parent_US_CNVs,
                                data.frame(SPARK_parent_US_SNVsExonicSizes_curr[[3]]))
}
SPARK_parent_US_CNVs$freq_max <- Determine_CNV_MAF(SPARK_parent_US_CNVs)
  
SPARK_parent_US_SNVs <- SPARK_parent_US_SNVs[SPARK_parent_US_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                       | SPARK_parent_US_SNVs$LoF,] 
SPARK_parent_US_SNVs$UID <- paste(SPARK_parent_US_SNVs$X.Sample, SPARK_parent_US_SNVs$X.id, sep='.')
SPARK_parent_US_SNVs <- SPARK_parent_US_SNVs[!(SPARK_parent_US_SNVs$LoF & SPARK_parent_US_SNVs$freq_max > 0.01), ]
SPARK_parent_US_SNVs <- SPARK_parent_US_SNVs[!is.na(SPARK_parent_US_SNVs$X.Sample),]

# SSC_parent_US_SV_SNVsExonicSizes <- yaml::yaml.load_file("./SV/SSC_parent_US_SV_SNVsExonicSizes.yaml")
# SSC_parent_US_SV_SNVs <- data.frame(SSC_parent_US_SV_SNVsExonicSizes[[1]])
# SSC_parent_US_SV_SNVs <- SSC_parent_US_SV_SNVs[SSC_parent_US_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
#                                                | SSC_parent_US_SV_SNVs$LoF,] 
# SSC_parent_US_SV_SNVs <- SSC_parent_US_SV_SNVs[SSC_parent_US_SV_SNVs$cnvENDAnn - SSC_parent_US_SV_SNVs$cnvSTARTAnn <= 10000, ]
# SSC_parent_US_SV_SNVs$UID <- paste(SSC_parent_US_SV_SNVs$X.Sample, SSC_parent_US_SV_SNVs$X.id, sep='.')
# SSC_parent_US_SV_SNVs <- SSC_parent_US_SV_SNVs[!(SSC_parent_US_SV_SNVs$LoF & SSC_parent_US_SV_SNVs$freq_max > 0.01), ]
# SSC_parent_US_SV_ExonicSizes <- data.frame(SSC_parent_US_SV_SNVsExonicSizes[[2]])
# SSC_parent_US_all_ExonicSizes <- merge(SSC_parent_US_ExonicSizes, SSC_parent_US_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T)
# #Combine
# SSC_parent_US_all_ExonicSizes <- merge(SSC_parent_US_ExonicSizes, SSC_parent_US_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
# SSC_parent_US_all_ExonicSizes$exonicSize.x[is.na(SSC_parent_US_all_ExonicSizes$exonicSize.x)] <- 0
# SSC_parent_US_all_ExonicSizes$exonicSize.y[is.na(SSC_parent_US_all_ExonicSizes$exonicSize.y)] <- 0
# SSC_parent_US_all_ExonicSizes$exonicSize <- pmax(SSC_parent_US_all_ExonicSizes$exonicSize.x, SSC_parent_US_all_ExonicSizes$exonicSize.y)
# SSC_parent_US_combined_SNVs <- rbind(SSC_parent_US_SNVs, SSC_parent_US_SV_SNVs)
# SSC_parent_US_proc_SNVs <- dplyr::distinct(SSC_parent_US_combined_SNVs, UID, .keep_all=T)
# 
# 
# # Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
# SSC_parent_US_SVs <- data.frame(SSC_parent_US_SV_SNVsExonicSizes[[3]])
# SSC_parent_US_SVs <- SSC_parent_US_SVs[SSC_parent_US_SVs$length <= 10000
#                                        & SSC_parent_US_SVs$Sample.ID %in% SSC_US_meta$Sample.ID, ]
# 
# SSC_US_homozyg_SVs <- SSC_parent_US_SVs[SSC_parent_US_SVs$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
# SSC_US_homozyg_SVs$SVUID <- with(SSC_US_homozyg_SVs, paste0(sample,CHROM,START,END))
# 
# SSC_parent_US_proc_SNVs <- SSC_parent_US_proc_SNVs[!with(SSC_parent_US_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% SSC_US_homozyg_SVs$SVUID,]

SPARK_parent_US_proc_SNVs <- SPARK_parent_US_SNVs
SPARK_parent_US_proc_SNVs$Min_Dist <- by(SPARK_parent_US_proc_SNVs, seq_len(nrow(SPARK_parent_US_proc_SNVs)),
                                       function(r) r$MIN = min(abs(r$POS - SPARK_parent_US_proc_SNVs$POS[SPARK_parent_US_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                         SPARK_parent_US_proc_SNVs$CHROM == r$CHROM &
                                                                                                         SPARK_parent_US_proc_SNVs$POS != r$POS])) > 50)
SPARK_parent_US_proc_SNVs <- SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$Min_Dist,]

SPARK_parent_US_clogit_res <- Get_Log_Reg(SPARK_US_meta,
                                        SPARK_parent_US_proc_SNVs,
                                        gene_set_data_path,
                                        SPARK_parent_US_ExonicSizes,
                                        ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                        mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                        type='clogit')
### SPARK Proband-Unaffected Siblings

SPARK_sibling_meta <- data.frame()
for (i in 1:3) {
  SPARK_sibling_metadata_path <- paste("./SPARK/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  SPARK_sibling_meta <- rbind(SPARK_sibling_meta, Get_Filtered_Metadata(SPARK_sibling_metadata_path, 
                                                                        child_relation = 'proband', 
                                                                        comparison='child-child'))
}

SPARK_proband_US_proc_SNVs <- rbind(SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$X.Sample %in% SPARK_sibling_meta$Sample.ID, ],
                                  SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$X.Sample %in% SPARK_sibling_meta$Sample.ID, ])
SPARK_proband_US_ExonicSizes <- rbind(SPARK_parent_US_ExonicSizes[SPARK_parent_US_ExonicSizes$Sample.ID %in% SPARK_sibling_meta$Sample.ID, ],
                                      SPARK_parent_proband_ExonicSizes[SPARK_parent_proband_ExonicSizes$Sample.ID %in% SPARK_sibling_meta$Sample.ID, ])

SPARK_proband_US_clogit_res <- Get_Log_Reg(SPARK_sibling_meta,
                                           SPARK_proband_US_proc_SNVs,
                                         gene_set_data_path,
                                         SPARK_proband_US_ExonicSizes,
                                         snv_count_groups = 'proband-unaffected sibling',
                                         ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                         mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                         type='clogit')

#### CH Event Limited to 0.1% *NOTE: Overwrites previous vars*####
### SPARK Parent-Proband
SPARK_parent_proband_0.1CH_SNVs <- SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$CHfreq <= 0.001, ]
SPARK_parent_proband_0.1CH_clogit_res <- Get_Log_Reg(SPARK_meta,
                                                   SPARK_parent_proband_0.1CH_SNVs,
                                                   gene_set_data_path,
                                                   SPARK_parent_proband_ExonicSizes,
                                                   ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                   mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                   type='clogit')

SPARK_parent_proband_0.01CH_SNVs <- SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$CHfreq <= 0.0001, ]
SPARK_parent_proband_0.01CH_clogit_res <- Get_Log_Reg(SPARK_meta,
                                                    SPARK_parent_proband_0.01CH_SNVs,
                                                    gene_set_data_path,
                                                    SPARK_parent_proband_ExonicSizes,
                                                    ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                    mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                    type='clogit')
### SPARK Parent-Unaffected Siblings
SPARK_parent_US_0.1CH_SNVs <- SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$CHfreq <= 0.001, ]
SPARK_parent_US_0.1CH_clogit_res <- Get_Log_Reg(SPARK_US_meta,
                                              SPARK_parent_US_0.1CH_SNVs,
                                              gene_set_data_path,
                                              SPARK_parent_US_ExonicSizes,
                                              ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                              mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                              type='clogit')

SPARK_parent_US_0.01CH_SNVs <- SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$CHfreq <= 0.0001, ]
SPARK_parent_US_0.01CH_clogit_res <- Get_Log_Reg(SPARK_US_meta,
                                               SPARK_parent_US_0.01CH_SNVs,
                                               gene_set_data_path,
                                               SPARK_parent_US_ExonicSizes,
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

### SPARK Proband-Unaffected Siblings
SPARK_proband_US_0.1CH_SNVs <- SPARK_proband_US_proc_SNVs[SPARK_proband_US_proc_SNVs$CHfreq <= 0.001, ]
SPARK_proband_US_0.1CH_clogit_res <- Get_Log_Reg(SPARK_sibling_meta,
                                               SPARK_proband_US_0.1CH_SNVs,
                                               gene_set_data_path,
                                               SPARK_proband_US_ExonicSizes,
                                               snv_count_groups = 'proband-unaffected sibling',
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

SPARK_proband_US_0.01CH_SNVs <- SPARK_proband_US_proc_SNVs[SPARK_proband_US_proc_SNVs$CHfreq <= 0.0001, ]
SPARK_proband_US_0.01CH_clogit_res <- Get_Log_Reg(SPARK_sibling_meta,
                                                SPARK_proband_US_0.01CH_SNVs,
                                                gene_set_data_path,
                                                SPARK_proband_US_ExonicSizes,
                                                snv_count_groups = 'proband-unaffected sibling',
                                                ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                type='clogit')

#### Different SNV MAF Tests ####
SPARK_parent_proband_proc_SNVs_10p <- SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$freq_max <= 0.1,]
SPARK_parent_proband_SNV10p_clogit_res <- Get_Log_Reg(SPARK_meta,
                                                    SPARK_parent_proband_proc_SNVs_10p,
                                                    gene_set_data_path,
                                                    SPARK_parent_proband_ExonicSizes,
                                                    ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                    mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                    type='clogit')
SPARK_parent_proband_proc_SNVs_1p <- SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$freq_max <= 0.01,]
SPARK_parent_proband_SNV1p_clogit_res <- Get_Log_Reg(SPARK_meta,
                                                   SPARK_parent_proband_proc_SNVs_1p,
                                                   gene_set_data_path,
                                                   SPARK_parent_proband_ExonicSizes,
                                                   ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                   mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                   type='clogit')

SPARK_parent_US_proc_SNVs_10p <- SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$freq_max <= 0.1,]
SPARK_parent_US_SNV10p_clogit_res <- Get_Log_Reg(SPARK_US_meta,
                                               SPARK_parent_US_proc_SNVs_10p,
                                               gene_set_data_path,
                                               SPARK_parent_US_ExonicSizes,
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')
SPARK_parent_US_proc_SNVs_1p <- SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$freq_max <= 0.01,]
SPARK_parent_US_SNV1p_clogit_res <- Get_Log_Reg(SPARK_US_meta,
                                              SPARK_parent_US_proc_SNVs_1p,
                                              gene_set_data_path,
                                              SPARK_parent_US_ExonicSizes,
                                              ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                              mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                              type='clogit')

SPARK_proband_US_proc_SNVs_10p <- SPARK_proband_US_proc_SNVs[SPARK_proband_US_proc_SNVs$freq_max <= 0.1,]
SPARK_proband_US_SNV10p_clogit_res <- Get_Log_Reg(SPARK_sibling_meta,
                                                SPARK_proband_US_proc_SNVs_10p,
                                                gene_set_data_path,
                                                SPARK_proband_US_ExonicSizes,
                                                snv_count_groups = 'proband-unaffected sibling',
                                                ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                type='clogit')
SPARK_proband_US_proc_SNVs_1p <- SPARK_proband_US_proc_SNVs[SPARK_proband_US_proc_SNVs$freq_max <= 0.01,]
SPARK_proband_US_SNV1p_clogit_res <- Get_Log_Reg(SPARK_sibling_meta,
                                               SPARK_proband_US_proc_SNVs_1p,
                                               gene_set_data_path,
                                               SPARK_proband_US_ExonicSizes,
                                               snv_count_groups = 'proband-unaffected sibling',
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')
#### Extra results and data formatting ####
Get_Regression_Slide_Data <- function(regression_results,save_path='../DT/LogRegResults/',
                                      include_base_reg=F) {
  res <- list()
  for (i in 1:length(regression_results)) {
    message(i)
    slide_data <- data.frame()
    result <- regression_results[[i]]
    name <- names(regression_results)[[i]]
    all_coeff <- result[[1]]$all[[2]]$coefficients[[length(result[[1]]$all[[2]]$coefficients)]]
    all_P <- result[[1]]$all[[3]][[length(result[[1]]$all[[3]])]][[2]] 
    all_G1 <- result[[1]]$all$Group1
    all_G2 <- result[[1]]$all$Group2
    all_G3 <- result[[1]]$all$Group3
    if (length(result[[1]]$syn) > 0) {
      syn_coeff <- result[[1]]$syn[[2]]$coefficients[[length(result[[1]]$syn[[2]]$coefficients)]]
      syn_P <- result[[1]]$syn[[3]][[length(result[[1]]$syn[[3]])]][[2]]
      syn_G1 <- result[[1]]$syn$Group1
      syn_G2 <- result[[1]]$syn$Group2
      syn_G3 <- result[[1]]$syn$Group3
    }
    else {
      syn_coeff <- NA
      syn_P <- NA
      syn_G1 <- NA
      syn_G2 <- NA
      syn_G3 <- NA
    }
    mis_coeff <- result[[1]]$mis[[2]]$coefficients[[length(result[[1]]$mis[[2]]$coefficients)]]
    mis_P <- result[[1]]$mis[[3]][[length(result[[1]]$mis[[3]])]][[2]]
    mis_G1 <- result[[1]]$mis$Group1
    mis_G2 <- result[[1]]$mis$Group2
    mis_G3 <- result[[1]]$mis$Group3
    if (length(result[[1]]$lof) > 0) {
      lof_coeff <- result[[1]]$lof[[2]]$coefficients[[length(result[[1]]$lof[[2]]$coefficients)]]
      lof_P <- result[[1]]$lof[[3]][[length(result[[1]]$lof[[3]])]][[2]]
      lof_G1 <- result[[1]]$lof$Group1
      lof_G2 <- result[[1]]$lof$Group2
      lof_G3 <- result[[1]]$lof$Group3
    }
    else {
      lof_coeff <- NA
      lof_P <- NA
      lof_G1 <- NA
      lof_G2 <- NA
      lof_G3 <- NA
    }
    neuro_coeff <- result[[2]]$Neurof_UnionInclusive[[2]]$coefficients[[length(result[[2]]$Neurof_UnionInclusive[[2]]$coefficients)]]
    neuro_P <- result[[2]]$Neurof_UnionInclusive[[3]][[length(result[[2]]$Neurof_UnionInclusive[[3]])]][[2]]
    neuro_G1 <- result[[2]]$Neurof_UnionInclusive$Group1
    neuro_G2 <- result[[2]]$Neurof_UnionInclusive$Group2
    neuro_G3 <- result[[2]]$Neurof_UnionInclusive$Group3
    # if (i == 2) {
    #   neuro_coeff <- 0
    #   neuro_P <- 0
    #   neuro_G1 <- 0
    #   neuro_G2 <- 0
    #   neuro_G3 <- 0
    # }
    if (!is.numeric(neuro_coeff)) {
      neuro_coeff <- NA
      neuro_P <- NA
      neuro_G1 <- NA
      neuro_G2 <- NA
      neuro_G3 <- NA
    }
    slide_data <- rbind(slide_data, 
                        data.frame(coeff=round(all_coeff,3), P=round(all_P,3), G1=all_G1, G2=all_G2, G3=all_G3),
                        data.frame(coeff=round(syn_coeff,3), P=round(syn_P,3), G1=syn_G1, G2=syn_G2, G3=syn_G3),
                        data.frame(coeff=round(mis_coeff,3), P=round(mis_P,3), G1=mis_G1, G2=mis_G2, G3=mis_G3),
                        data.frame(coeff=round(lof_coeff,3), P=round(lof_P,3), G1=lof_G1, G2=lof_G2, G3=lof_G3),
                        data.frame(coeff=round(neuro_coeff,3), P=round(neuro_P,3), G1=neuro_G1, G2=neuro_G2, G3=neuro_G3))
    if (include_base_reg) {
      all_base_coeff <- result[[1]]$all_base[[2]]$coefficients[[length(result[[1]]$all_base[[2]]$coefficients)]]
      all_base_P <- result[[1]]$all_base[[3]][[length(result[[1]]$all_base[[3]])]][[2]]
      all_base_G1 <- result[[1]]$all_base$Group1
      all_base_G2 <- result[[1]]$all_base$Group2
      all_base_G3 <- result[[1]]$all_base$Group3
      slide_data <- rbind(data.frame(coeff=round(all_base_coeff,3), P=round(all_base_P,3), G1=all_base_G1, G2=all_base_G2, G3=all_base_G3), 
                          slide_data)
    }
    
    write.csv(slide_data, paste(paste(save_path,name,sep=''),'.csv',sep=''))
    res[[name]] <- slide_data
  }
  
  return (res)
}

SPARK_reg_resuts <- list(SPARK_parent_proband_clogit_res=SPARK_parent_proband_clogit_res,
                   SPARK_parent_US_clogit_res=SPARK_parent_US_clogit_res,
                   SPARK_proband_US_clogit_res=SPARK_proband_US_clogit_res,
                   SPARK_parent_proband_0.1CH_clogit_res=SPARK_parent_proband_0.1CH_clogit_res,
                   SPARK_parent_US_0.1CH_clogit_res=SPARK_parent_US_0.1CH_clogit_res,
                   SPARK_proband_US_0.1CH_clogit_res=SPARK_proband_US_0.1CH_clogit_res,
                   SPARK_parent_proband_0.01CH_clogit_res=SPARK_parent_proband_0.01CH_clogit_res,
                   SPARK_parent_US_0.01CH_clogit_res=SPARK_parent_US_0.01CH_clogit_res,
                   SPARK_proband_US_0.01CH_clogit_res=SPARK_proband_US_0.01CH_clogit_res,
                   SPARK_parent_proband_SNV10p_clogit_res=SPARK_parent_proband_SNV10p_clogit_res,
                   SPARK_parent_proband_SNV1p_clogit_res=SPARK_parent_proband_SNV1p_clogit_res,
                   SPARK_parent_US_SNV10p_clogit_res=SPARK_parent_US_SNV10p_clogit_res,
                   SPARK_parent_US_SNV1p_clogit_res=SPARK_parent_US_SNV1p_clogit_res,
                   SPARK_proband_US_SNV10p_clogit_res=SPARK_proband_US_SNV10p_clogit_res,
                   SPARK_proband_US_SNV1p_clogit_res=SPARK_proband_US_SNV1p_clogit_res)

SPARK_slide_data <- format(Get_Regression_Slide_Data(SPARK_reg_resuts, include_base_reg=F), 
                     scientific=F)


#SNV Counts for all SV:
SPARK_all_SNVs_combined <- rbind(SPARK_parent_proband_proc_SNVs, SPARK_parent_US_proc_SNVs)
SPARK_all_SNVs_combined_unique <- dplyr::distinct(SPARK_all_SNVs_combined, UID, .keep_all = T)
SPARK_all_ExonicSizes_combined <- rbind(SPARK_parent_proband_ExonicSizes, SPARK_parent_US_ExonicSizes)
SPARK_all_ExonicSizes_combined_unique <- dplyr::distinct(SPARK_all_ExonicSizes_combined, Sample.ID, .keep_all = T)

#SNV Counts for all CNV
SNV_in_CNV_only_SPARK <- SPARK_all_SNVs_combined_unique[SPARK_all_SNVs_combined_unique$UID %in% SPARK_parent_proband_SNVs$UID |
                                                      SPARK_all_SNVs_combined_unique$UID %in% SPARK_parent_US_SNVs$UID, ]
all_CNV_count_SPARK <- length(SNV_in_CNV_only_SPARK$X.Sample)
syn_CNV_count_SPARK <- length(SNV_in_CNV_only_SPARK$X.Sample[SNV_in_CNV_only_SPARK$effect_priority == 'synonymous SNV'])
mis_CNV_count_SPARK <- length(SNV_in_CNV_only_SPARK$X.Sample[SNV_in_CNV_only_SPARK$effect_priority == 'nonsynonymous SNV'])
lof_CNV_count_SPARK <- length(SNV_in_CNV_only_SPARK$X.Sample[which(SNV_in_CNV_only_SPARK$LoF)])

# Other SNV in CNV/SV count metrics
avg_bp_cnv_SPARK <- summary(SPARK_all_ExonicSizes_combined_unique$exonicSize[SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SPARK$X.Sample])
snv_per_bp_CNV_SPARK <- all_CNV_count_SPARK/ sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SPARK$X.Sample])

# No. of CNV & SV counts
SPARK_parent_proband_CNVs <- SPARK_parent_proband_CNVs[SPARK_parent_proband_CNVs$freq_max <= 0.01
                                                   & SPARK_parent_proband_CNVs$Sample.ID %in% SPARK_meta$Sample.ID, ]
SPARK_parent_US_CNVs <- SPARK_parent_US_CNVs[SPARK_parent_US_CNVs$freq_max <= 0.01
                                         & SPARK_parent_US_CNVs$Sample.ID %in% SPARK_US_meta$Sample.ID, ]
SPARK_parent_child_CNVs <- rbind(SPARK_parent_proband_CNVs, SPARK_parent_US_CNVs)
SPARK_parent_child_CNVs <- SPARK_parent_child_CNVs[!duplicated(paste(SPARK_parent_child_CNVs$Sample.ID, 
                                                                 SPARK_parent_child_CNVs$CHROM,
                                                                 SPARK_parent_child_CNVs$START,
                                                                 SPARK_parent_child_CNVs$END, sep='.')), ]

SPARK_father_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'father']
SPARK_mother_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'mother']
SPARK_proband_IDs <- SPARK_meta$Sample.ID[SPARK_meta$Relation == 'proband']
SPARK_father_US_IDs <- SPARK_US_meta$Sample.ID[SPARK_US_meta$Relation == 'father']
SPARK_mother_US_IDs <- SPARK_US_meta$Sample.ID[SPARK_US_meta$Relation == 'mother']
SPARK_US_IDs <- SPARK_US_meta$Sample.ID[!SPARK_US_meta$Relation %in% c('father', 'mother')]

SPARK_CNV_count <- nrow(SPARK_parent_child_CNVs)
SPARK_CNV_father <- SPARK_parent_proband_CNVs[SPARK_parent_proband_CNVs$Sample.ID %in% SPARK_father_IDs, ]
SPARK_CNV_count_father <- nrow(SPARK_CNV_father)
SPARK_CNV_mother <- SPARK_parent_proband_CNVs[SPARK_parent_proband_CNVs$Sample.ID %in% SPARK_mother_IDs, ]
SPARK_CNV_count_mother <- nrow(SPARK_CNV_mother)
SPARK_CNV_proband <- SPARK_parent_proband_CNVs[SPARK_parent_proband_CNVs$Sample.ID %in% SPARK_proband_IDs, ]
SPARK_CNV_count_proband <- nrow(SPARK_CNV_proband)
SPARK_CNV_father_US <- SPARK_parent_US_CNVs[SPARK_parent_US_CNVs$Sample.ID %in% SPARK_father_US_IDs, ]
SPARK_CNV_count_father_US <- nrow(SPARK_CNV_father_US)
SPARK_CNV_mother_US <- SPARK_parent_US_CNVs[SPARK_parent_US_CNVs$Sample.ID %in% SPARK_mother_US_IDs, ]
SPARK_CNV_count_mother_US <- nrow(SPARK_CNV_mother_US)
SPARK_CNV_US <- SPARK_parent_US_CNVs[SPARK_parent_US_CNVs$Sample.ID %in% SPARK_US_IDs, ]
SPARK_CNV_count_US <- nrow(SPARK_CNV_US)

## Total exon lengths from each individual group

SPARK_CNV_father_totalExonLen <- sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[
  SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SPARK_CNV_father$Sample.ID])
SPARK_CNV_mother_totalExonLen <- sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[
  SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SPARK_CNV_mother$Sample.ID])
SPARK_CNV_proband_totalExonLen <- sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[
  SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SPARK_CNV_proband$Sample.ID])
SPARK_CNV_father_US_totalExonLen <- sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[
  SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SPARK_CNV_father_US$Sample.ID])
SPARK_CNV_mother_US_totalExonLen <- sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[
  SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SPARK_CNV_mother_US$Sample.ID])
SPARK_CNV_US_totalExonLen <- sum(SPARK_all_ExonicSizes_combined_unique$exonicSize[
  SPARK_all_ExonicSizes_combined_unique$Sample.ID %in% SPARK_CNV_US$Sample.ID])

Get_Number_Of_Individual_Comparison <- function(data, meta,
                                                comparison='parent-child',
                                                child_relation='proband') {
  res <- data.frame()
  total_res <- data.frame()
  child_relations <- ifelse(child_relation=='proband', 
                            c('proband'), c('unaffected sibling', 'other sibling'))
  child_rel_label <- ifelse(child_relation=='proband', 'Probands', 'Unaffected Siblings')
  total=length(unique(meta$Sample.ID[meta$Relation %in% child_relations],))
  all_SNVs_child <- data[data$X.Sample %in% meta$Sample.ID[meta$Relation %in% child_relations], ]
  all_SNVs_child_count <- length(unique(all_SNVs_child$X.Sample)) / total
  syn_SNVs_child_count <- length(unique(all_SNVs_child$X.Sample[all_SNVs_child$effect_priority == 'synonymous SNV'])) / total
  mis_SNVs_child_count <- length(unique(all_SNVs_child$X.Sample[all_SNVs_child$effect_priority == 'nonsynonymous SNV'])) / total
  lof_SNVs_child_count <- length(unique(all_SNVs_child$X.Sample[which(all_SNVs_child$LoF)])) / total
  
  res <- rbind(res, data.frame(Relation=child_rel_label,
                               `Variant Type`='All Variants',
                               Count=all_SNVs_child_count),
               data.frame(Relation=child_rel_label,
                          `Variant Type`='Synonymous',
                          Count=syn_SNVs_child_count),
               data.frame(Relation=child_rel_label,
                          `Variant Type`='Missense',
                          Count=mis_SNVs_child_count),
               data.frame(Relation=child_rel_label,
                          `Variant Type`='LoF',
                          Count=lof_SNVs_child_count))
  total_res <- rbind(total_res, data.frame(Relation=child_rel_label,
                                           Count=total))
  if (comparison=='parent-child') {
    total=length(unique(meta$Sample.ID[meta$Relation=='father']))
    all_SNVs_father <- data[data$X.Sample %in% meta$Sample.ID[meta$Relation == 'father'], ]
    all_SNVs_father_count <- length(unique(all_SNVs_father$X.Sample)) / total
    syn_SNVs_father_count <- length(unique(all_SNVs_father$X.Sample[all_SNVs_father$effect_priority == 'synonymous SNV'])) / total
    mis_SNVs_father_count <- length(unique(all_SNVs_father$X.Sample[all_SNVs_father$effect_priority == 'nonsynonymous SNV'])) / total
    lof_SNVs_father_count <- length(unique(all_SNVs_father$X.Sample[which(all_SNVs_father$LoF)])) / total
    
    res <- rbind(res, data.frame(Relation='Fathers',
                                 `Variant Type`='All Variants',
                                 Count=all_SNVs_father_count),
                 data.frame(Relation='Fathers',
                            `Variant Type`='Synonymous',
                            Count=syn_SNVs_father_count),
                 data.frame(Relation='Fathers',
                            `Variant Type`='Missense',
                            Count=mis_SNVs_father_count),
                 data.frame(Relation='Fathers',
                            `Variant Type`='LoF',
                            Count=lof_SNVs_father_count))
    total_res <- rbind(total_res, data.frame(Relation='Fathers',
                                             Count=total))
    total=length(unique(meta$Sample.ID[meta$Relation=='mother']))
    all_SNVs_mother <- data[data$X.Sample %in% meta$Sample.ID[meta$Relation == 'mother'], ]
    all_SNVs_mother_count <- length(unique(all_SNVs_mother$X.Sample)) / total
    syn_SNVs_mother_count <- length(unique(all_SNVs_mother$X.Sample[all_SNVs_mother$effect_priority == 'synonymous SNV'])) / total
    mis_SNVs_mother_count <- length(unique(all_SNVs_mother$X.Sample[all_SNVs_mother$effect_priority == 'nonsynonymous SNV'])) / total
    lof_SNVs_mother_count <- length(unique(all_SNVs_mother$X.Sample[which(all_SNVs_mother$LoF)])) / total
    
    res <- rbind(res, data.frame(Relation='Mothers',
                                 `Variant Type`='All Variants',
                                 Count=all_SNVs_mother_count),
                 data.frame(Relation='Mothers',
                            `Variant Type`='Synonymous',
                            Count=syn_SNVs_mother_count),
                 data.frame(Relation='Mothers',
                            `Variant Type`='Missense',
                            Count=mis_SNVs_mother_count),
                 data.frame(Relation='Mothers',
                            `Variant Type`='LoF',
                            Count=lof_SNVs_mother_count))
    total_res <- rbind(total_res, data.frame(Relation='Mothers',
                                             Count=total))
  }
  else if (comparison == 'child-child') {
    total=length(unique(meta$Sample.ID[meta$Relation%in% c('unaffected sibling', 'other sibling')]))
    all_SNVs_US <- data[data$X.Sample %in% meta$Sample.ID[meta$Relation %in% c('unaffected sibling', 'other sibling')], ]
    all_SNVs_US_count <- length(unique(all_SNVs_US$X.Sample)) / total
    syn_SNVs_US_count <- length(unique(all_SNVs_US$X.Sample[all_SNVs_US$effect_priority == 'synonymous SNV'])) / total
    mis_SNVs_US_count <- length(unique(all_SNVs_US$X.Sample[all_SNVs_US$effect_priority == 'nonsynonymous SNV'])) / total
    lof_SNVs_US_count <- length(unique(all_SNVs_US$X.Sample[which(all_SNVs_US$LoF)])) / total
    
    res <- rbind(res, data.frame(Relation='Unaffected Siblings',
                                 `Variant Type`='All Variants',
                                 Count=all_SNVs_US_count),
                 data.frame(Relation='Unaffected Siblings',
                            `Variant Type`='Synonymous',
                            Count=syn_SNVs_US_count),
                 data.frame(Relation='Unaffected Siblings',
                            `Variant Type`='Missense',
                            Count=mis_SNVs_US_count),
                 data.frame(Relation='Unaffected Siblings',
                            `Variant Type`='LoF',
                            Count=lof_SNVs_US_count))
    total_res <- rbind(total_res, data.frame(Relation='Unaffected Siblings',
                                             Count=total))
  }
  
  return (list(res, total_res))
}

Generate_Indiv_Count_Bar_Plots <- function(data, name, save_path='../DT/IndivCountBarPlots/',
                                           comparison_type='parent-child-proband') {
  if (comparison_type == 'parent-child-proband') {
    legend_level = c('Probands', 'Fathers', 'Mothers')
  } else if (comparison_type == 'parent-child-US') {
    legend_level = c('Unaffected Siblings', 'Fathers', 'Mothers')
  } else if (comparison_type == 'child-child') {
    legend_level = c('Probands', 'Unaffected Siblings')
  }
  ggplot(data, 
         aes(x=factor(Variant.Type, level=c('All Variants', 'Synonymous', 'Missense', 'LoF')),
             y=Count,fill=factor(Relation, level=legend_level))) + 
    geom_bar(stat='identity',position='dodge') + 
    geom_text(aes(label=round(Count,3)), position=position_dodge(width=0.9), vjust=-0.25) + 
    xlab('Variant Type') + ylab('Normalized Count') + labs(fill='Relation') + 
    scale_fill_hue(l=45)
  
  ggsave(paste(save_path, paste(name, '.png',sep=''), sep=''))
}

SPARK_parent_proband_indiv_data <- Get_Number_Of_Individual_Comparison(SPARK_parent_proband_proc_SNVs, 
                                                                     SPARK_meta,
                                                                     comparison='parent-child', 
                                                                     child_relation='proband')
SPARK_parent_US_indiv_data <- Get_Number_Of_Individual_Comparison(SPARK_parent_US_proc_SNVs, 
                                                                SPARK_US_meta,
                                                                comparison='parent-child', 
                                                                child_relation='unaffected sibling')
SPARK_proband_US_indiv_data <- Get_Number_Of_Individual_Comparison(SPARK_proband_US_proc_SNVs, 
                                                                 SPARK_sibling_meta,
                                                                 comparison='child-child', 
                                                                 child_relation='proband')

SPARK_parent_proband_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(SPARK_parent_proband_0.1CH_SNVs, 
                                                                         SPARK_meta,
                                                                         comparison='parent-child', 
                                                                         child_relation='proband')
SPARK_parent_US_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(SPARK_parent_US_0.1CH_SNVs, 
                                                                    SPARK_US_meta,
                                                                    comparison='parent-child', 
                                                                    child_relation='unaffected sibling')
SPARK_proband_US_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(SPARK_proband_US_0.1CH_SNVs, 
                                                                     SPARK_sibling_meta,
                                                                     comparison='child-child', 
                                                                     child_relation='proband')

Generate_Indiv_Count_Bar_Plots(SPARK_parent_proband_indiv_data[[1]], 'SPARK_parent_proband_indiv_data')
Generate_Indiv_Count_Bar_Plots(SPARK_parent_US_indiv_data[[1]], 'SPARK_parent_US_indiv_data', comparison_type='parent-child-US')
Generate_Indiv_Count_Bar_Plots(SPARK_proband_US_indiv_data[[1]], 'SPARK_proband_US_indiv_data', comparison_type='child-child')
Generate_Indiv_Count_Bar_Plots(SPARK_parent_proband_indiv_0.1_data[[1]], 'SPARK_parent_proband_indiv_0.1_data')
Generate_Indiv_Count_Bar_Plots(SPARK_parent_US_indiv_0.1_data[[1]], 'SPARK_parent_US_indiv_0.1_data', comparison_type='parent-child-US')
Generate_Indiv_Count_Bar_Plots(SPARK_proband_US_indiv_0.1_data[[1]], 'SPARK_proband_US_indiv_0.1_data', comparison_type='child-child')

## Create list of all gene set regression results
Get_Gene_Set_Regression_Data <- function(regression_results, P_cutoff=0.1,
                                         save_path='../DT/GSLogRegResults/') {
  res <- list()
  exclude <- c('GOP_lessthan500g', 'PhMm_CompleteLethality')
  for (i in 1:length(regression_results)) {
    message(i)
    slide_data <- data.frame()
    result <- regression_results[[i]]
    name <- names(regression_results)[[i]]
    gs <- result[[2]]
    for (j in 1:length(gs)) {
      curr_gs <- gs[[j]]
      gs_name <- names(gs)[[j]]
      coeff <- curr_gs[[2]]$coefficients[[length(curr_gs[[2]]$coefficients)]]
      P <- curr_gs[[3]][[length(curr_gs[[3]])]][[2]]
      if (!gs_name %in% exclude & !grepl('scRNA',gs_name,fixed=T) & P <= P_cutoff) {
        G1 <- curr_gs$Group1
        G2 <- curr_gs$Group2
        G3 <- curr_gs$Group3
        slide_data <- rbind(slide_data, 
                            data.frame(name=gs_name,coeff=round(coeff,3), P=round(P,3), G1=G1, G2=G2, G3=G3))
        
      }
    }
    slide_data <- slide_data[order(slide_data$P),]
    write.csv(slide_data, paste(paste(save_path,name,sep=''),'_genesets.csv',sep=''))
    res[[name]] <- slide_data
  }
  
  return (res)
}
SPARK_gs_reg_result_list <- list(SPARK_parent_proband_0.1CH_clogit_res=SPARK_parent_proband_0.1CH_clogit_res,
                           SPARK_parent_US_0.1CH_clogit_res=SPARK_parent_US_0.1CH_clogit_res,
                           SPARK_proband_US_0.1CH_clogit_res=SPARK_proband_US_0.1CH_clogit_res)

SPARK_gs_reg_data <- format(Get_Gene_Set_Regression_Data(SPARK_gs_reg_result_list), scientific=F)

## Genes with multiple hits/multiple samples
SPARK_gene_summary <- SPARK_all_SNVs_combined_unique %>% filter(X.Sample%in%SPARK_proband_IDs & effect_priority=='nonsynonymous SNV' & gnomAD_pRec >= 0.9 & gnomAD_oe_lof_upper >= 0.35) %>% 
  group_by(gene_symbol) %>% summarise(Proband=length(unique(SPARK_all_SNVs_combined_unique$X.Sample[SPARK_all_SNVs_combined_unique$X.Sample %in% SPARK_proband_IDs & SPARK_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & SPARK_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV'])), 
                                      Father=length(unique(SPARK_all_SNVs_combined_unique$X.Sample[SPARK_all_SNVs_combined_unique$X.Sample %in% SPARK_father_IDs & SPARK_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & SPARK_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV'])),
                                      Mother=length(unique(SPARK_all_SNVs_combined_unique$X.Sample[SPARK_all_SNVs_combined_unique$X.Sample %in% SPARK_mother_IDs & SPARK_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & SPARK_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV'])),
                                      US=length(unique(SPARK_all_SNVs_combined_unique$X.Sample[SPARK_all_SNVs_combined_unique$X.Sample %in% SPARK_US_IDs & SPARK_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & SPARK_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV']))) %>% arrange(-Proband)
SPARK_gene_summary <- data.frame(SPARK_gene_summary[SPARK_gene_summary$Proband > 1,])
write.csv(SPARK_gene_summary, '../DT/SPARK_proband_genes_missense_summary.csv')

## Gene set enrichment analysis
# MSSNG_proband_genes_mis <- MSSNG_parent_proband_proc_SNVs$gene_symbol[MSSNG_parent_proband_proc_SNVs$effect_priority == 'nonsynonymous SNV'
#                                                                       & MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_proband_IDs
#                                                                       & MSSNG_parent_proband_proc_SNVs$gnomAD_pRec >= 0.9 
#                                                                       & MSSNG_parent_proband_proc_SNVs$gnomAD_oe_lof_upper >= 0.35]
# write(unique(MSSNG_proband_genes_mis),
#       '../DT/MSSNG_proband_genes_missense.txt')
# MSSNG_proband_genes_mis_dt <- data.frame(table(MSSNG_proband_genes_mis))
# 
# SSC_proband_genes_mis <- SSC_parent_proband_proc_SNVs$gene_symbol[SSC_parent_proband_proc_SNVs$effect_priority == 'nonsynonymous SNV'
#                                                                   & SSC_parent_proband_proc_SNVs$X.Sample %in% SSC_proband_IDs
#                                                                   & SSC_parent_proband_proc_SNVs$gnomAD_pRec >= 0.9 
#                                                                   & SSC_parent_proband_proc_SNVs$gnomAD_oe_lof_upper >= 0.35]
# write(unique(SSC_proband_genes_mis),
#       '../DT/SSC_proband_genes_missense.txt')
# SSC_proband_genes_mis_dt <- data.frame(table(SSC_proband_genes_mis))

