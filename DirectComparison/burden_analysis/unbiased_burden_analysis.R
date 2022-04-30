library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")


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
### MSSNG Parent-Proband
MSSNG_metadata_path <- "MSSNG_metadata.tsv"
MSSNG_meta <- Get_Filtered_Metadata(MSSNG_metadata_path, child_relation = 'proband')
MSSNG_parent_proband_SNVsExonicSizes <- yaml::yaml.load_file("./SV/MSSNG_parent_proband_SNVsExonicSizes.yaml")
MSSNG_parent_proband_SNVs <- data.frame(MSSNG_parent_proband_SNVsExonicSizes[[1]])
MSSNG_parent_proband_SNVs <- MSSNG_parent_proband_SNVs[MSSNG_parent_proband_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                       | MSSNG_parent_proband_SNVs$LoF,] 
MSSNG_parent_proband_SNVs$UID <- paste(MSSNG_parent_proband_SNVs$X.Sample, MSSNG_parent_proband_SNVs$X.id, sep='.')
MSSNG_parent_proband_SNVs <- MSSNG_parent_proband_SNVs[!(MSSNG_parent_proband_SNVs$LoF & MSSNG_parent_proband_SNVs$freq_max > 0.01), ]
MSSNG_parent_proband_ExonicSizes <- data.frame(MSSNG_parent_proband_SNVsExonicSizes[[2]])
MSSNG_parent_proband_SV_SNVsExonicSizes <- yaml::yaml.load_file("./SV/MSSNG_parent_proband_SV_SNVsExonicSizes.yaml")
MSSNG_parent_proband_SV_SNVs <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[1]])
MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[MSSNG_parent_proband_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                             | MSSNG_parent_proband_SV_SNVs$LoF,] 
MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[MSSNG_parent_proband_SV_SNVs$cnvENDAnn - MSSNG_parent_proband_SV_SNVs$cnvSTARTAnn <= 10000, ]
MSSNG_parent_proband_SV_SNVs$UID <- paste(MSSNG_parent_proband_SV_SNVs$X.Sample, MSSNG_parent_proband_SV_SNVs$X.id, sep='.')
MSSNG_parent_proband_SV_SNVs <- MSSNG_parent_proband_SV_SNVs[!(MSSNG_parent_proband_SV_SNVs$LoF & MSSNG_parent_proband_SV_SNVs$freq_max > 0.01), ]
MSSNG_parent_proband_SV_ExonicSizes <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[2]])

#Combine
MSSNG_parent_proband_all_ExonicSizes <- merge(MSSNG_parent_proband_ExonicSizes, MSSNG_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
MSSNG_parent_proband_all_ExonicSizes$exonicSize.x[is.na(MSSNG_parent_proband_all_ExonicSizes$exonicSize.x)] <- 0
MSSNG_parent_proband_all_ExonicSizes$exonicSize.y[is.na(MSSNG_parent_proband_all_ExonicSizes$exonicSize.y)] <- 0
MSSNG_parent_proband_all_ExonicSizes$exonicSize <- pmax(MSSNG_parent_proband_all_ExonicSizes$exonicSize.x, MSSNG_parent_proband_all_ExonicSizes$exonicSize.y)
MSSNG_parent_proband_combined_SNVs <- rbind(MSSNG_parent_proband_SNVs, MSSNG_parent_proband_SV_SNVs)
MSSNG_parent_proband_proc_SNVs <- dplyr::distinct(MSSNG_parent_proband_combined_SNVs, UID, .keep_all=T)

# Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
MSSNG_parent_proband_SVs_ILMN <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[3]])
MSSNG_parent_proband_SVs_ILMN <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$length <= 10000
                                                               & MSSNG_parent_proband_SVs_ILMN$Sample.ID %in% MSSNG_meta$Sample.ID, ]
MSSNG_parent_proband_SVs_CG <- data.frame(MSSNG_parent_proband_SV_SNVsExonicSizes[[4]])
MSSNG_parent_proband_SVs_CG <- MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$length <= 10000
                                                           & MSSNG_parent_proband_SVs_CG$Sample.ID %in% MSSNG_meta$Sample.ID, ]
MSSNG_homozyg_SVs <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
MSSNG_homozyg_SVs <- rbind(MSSNG_homozyg_SVs,MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')])
MSSNG_homozyg_SVs$SVUID <- with(MSSNG_homozyg_SVs, paste0(sample,CHROM,START,END))

MSSNG_parent_proband_proc_SNVs <- MSSNG_parent_proband_proc_SNVs[!with(MSSNG_parent_proband_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% MSSNG_homozyg_SVs$SVUID,]

MSSNG_parent_proband_proc_SNVs$Min_Dist <- by(MSSNG_parent_proband_proc_SNVs, seq_len(nrow(MSSNG_parent_proband_proc_SNVs)),
   function(r) r$MIN = min(abs(r$POS - MSSNG_parent_proband_proc_SNVs$POS[MSSNG_parent_proband_proc_SNVs$X.Sample == r$X.Sample & 
                                                                            MSSNG_parent_proband_proc_SNVs$CHROM == r$CHROM &
                                                                            MSSNG_parent_proband_proc_SNVs$POS != r$POS])) > 50)
MSSNG_parent_proband_proc_SNVs <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$Min_Dist,]



MSSNG_parent_proband_clogit_res <- Get_Log_Reg(MSSNG_meta,
                                               MSSNG_parent_proband_proc_SNVs,
                                               gene_set_data_path,
                                               MSSNG_parent_proband_all_ExonicSizes,
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

# Testing only simplex families in MSSNG
MSSNG_meta_SPX <- MSSNG_meta[MSSNG_meta$Family.type == 'SPX',]
MSSNG_parent_proband_proc_SNVs_SPX <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_meta_SPX$Sample.ID,]
MSSNG_meta_MPX <- MSSNG_meta[MSSNG_meta$Family.type == 'MPX',]
MSSNG_parent_proband_proc_SNVs_MPX <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_meta_MPX$Sample.ID,]

### SSC Parent-Proband
SSC_metadata_path <- "SSC_metadata.tsv"
SSC_meta <- Get_Filtered_Metadata(SSC_metadata_path, child_relation = 'proband')

SSC_parent_proband_SNVsExonicSizes <- yaml::yaml.load_file("./SV/SSC_parent_proband_SNVsExonicSizes.yaml")
SSC_parent_proband_SNVs <- data.frame(SSC_parent_proband_SNVsExonicSizes[[1]])
SSC_parent_proband_SNVs <- SSC_parent_proband_SNVs[SSC_parent_proband_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                             | SSC_parent_proband_SNVs$LoF,] 
SSC_parent_proband_SNVs$UID <- paste(SSC_parent_proband_SNVs$X.Sample, SSC_parent_proband_SNVs$X.id, sep='.')
SSC_parent_proband_SNVs <- SSC_parent_proband_SNVs[!(SSC_parent_proband_SNVs$LoF & SSC_parent_proband_SNVs$freq_max > 0.01), ]
SSC_parent_proband_ExonicSizes <- data.frame(SSC_parent_proband_SNVsExonicSizes[[2]])
SSC_parent_proband_SV_SNVsExonicSizes <- yaml::yaml.load_file("./SV/SSC_parent_proband_SV_SNVsExonicSizes.yaml")
SSC_parent_proband_SV_SNVs <- data.frame(SSC_parent_proband_SV_SNVsExonicSizes[[1]])
SSC_parent_proband_SV_SNVs <- SSC_parent_proband_SV_SNVs[SSC_parent_proband_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                   | SSC_parent_proband_SV_SNVs$LoF,] 
SSC_parent_proband_SV_SNVs <- SSC_parent_proband_SV_SNVs[SSC_parent_proband_SV_SNVs$cnvENDAnn - SSC_parent_proband_SV_SNVs$cnvSTARTAnn <= 10000, ]
SSC_parent_proband_SV_SNVs$UID <- paste(SSC_parent_proband_SV_SNVs$X.Sample, SSC_parent_proband_SV_SNVs$X.id, sep='.')
SSC_parent_proband_SV_SNVs <- SSC_parent_proband_SV_SNVs[!(SSC_parent_proband_SV_SNVs$LoF & SSC_parent_proband_SV_SNVs$freq_max > 0.01), ]
SSC_parent_proband_SV_ExonicSizes <- data.frame(SSC_parent_proband_SV_SNVsExonicSizes[[2]])
SSC_parent_proband_all_ExonicSizes <- merge(SSC_parent_proband_ExonicSizes, SSC_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T)
#Combine
SSC_parent_proband_all_ExonicSizes <- merge(SSC_parent_proband_ExonicSizes, SSC_parent_proband_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
SSC_parent_proband_all_ExonicSizes$exonicSize.x[is.na(SSC_parent_proband_all_ExonicSizes$exonicSize.x)] <- 0
SSC_parent_proband_all_ExonicSizes$exonicSize.y[is.na(SSC_parent_proband_all_ExonicSizes$exonicSize.y)] <- 0
SSC_parent_proband_all_ExonicSizes$exonicSize <- pmax(SSC_parent_proband_all_ExonicSizes$exonicSize.x, SSC_parent_proband_all_ExonicSizes$exonicSize.y)
SSC_parent_proband_combined_SNVs <- rbind(SSC_parent_proband_SNVs, SSC_parent_proband_SV_SNVs)
SSC_parent_proband_proc_SNVs <- dplyr::distinct(SSC_parent_proband_combined_SNVs, UID, .keep_all=T)


# Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
SSC_parent_proband_SVs <- data.frame(SSC_parent_proband_SV_SNVsExonicSizes[[3]])
SSC_parent_proband_SVs <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$length <= 10000
                                                 & SSC_parent_proband_SVs$Sample.ID %in% SSC_meta$Sample.ID, ]

SSC_homozyg_SVs <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
SSC_homozyg_SVs$SVUID <- with(SSC_homozyg_SVs, paste0(sample,CHROM,START,END))

SSC_parent_proband_proc_SNVs <- SSC_parent_proband_proc_SNVs[!with(SSC_parent_proband_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% SSC_homozyg_SVs$SVUID,]

SSC_parent_proband_proc_SNVs$Min_Dist <- by(SSC_parent_proband_proc_SNVs, seq_len(nrow(SSC_parent_proband_proc_SNVs)),
                                              function(r) r$MIN = min(abs(r$POS - SSC_parent_proband_proc_SNVs$POS[SSC_parent_proband_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                                     SSC_parent_proband_proc_SNVs$CHROM == r$CHROM &
                                                                                                                     SSC_parent_proband_proc_SNVs$POS != r$POS])) > 50)
SSC_parent_proband_proc_SNVs <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$Min_Dist,]

SSC_parent_proband_clogit_res <- Get_Log_Reg(SSC_meta,
                                             SSC_parent_proband_proc_SNVs,
                                             gene_set_data_path,
                                             SSC_parent_proband_all_ExonicSizes,
                                             ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                             mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                             type='clogit')
### SSC Parent-Unaffected Siblings
SSC_US_metadata_path <- "SSC_metadata.tsv"
SSC_US_meta <- Get_Filtered_Metadata(SSC_US_metadata_path, child_relation = 'unaffected sibling')

SSC_parent_US_SNVsExonicSizes <- yaml::yaml.load_file("./SV/SSC_parent_US_SNVsExonicSizes.yaml")
SSC_parent_US_SNVs <- data.frame(SSC_parent_US_SNVsExonicSizes[[1]])
SSC_parent_US_SNVs <- SSC_parent_US_SNVs[SSC_parent_US_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                         | SSC_parent_US_SNVs$LoF,] 
SSC_parent_US_SNVs$UID <- paste(SSC_parent_US_SNVs$X.Sample, SSC_parent_US_SNVs$X.id, sep='.')
SSC_parent_US_SNVs <- SSC_parent_US_SNVs[!(SSC_parent_US_SNVs$LoF & SSC_parent_US_SNVs$freq_max > 0.01), ]
SSC_parent_US_ExonicSizes <- data.frame(SSC_parent_US_SNVsExonicSizes[[2]])
SSC_parent_US_SV_SNVsExonicSizes <- yaml::yaml.load_file("./SV/SSC_parent_US_SV_SNVsExonicSizes.yaml")
SSC_parent_US_SV_SNVs <- data.frame(SSC_parent_US_SV_SNVsExonicSizes[[1]])
SSC_parent_US_SV_SNVs <- SSC_parent_US_SV_SNVs[SSC_parent_US_SV_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                         | SSC_parent_US_SV_SNVs$LoF,] 
SSC_parent_US_SV_SNVs <- SSC_parent_US_SV_SNVs[SSC_parent_US_SV_SNVs$cnvENDAnn - SSC_parent_US_SV_SNVs$cnvSTARTAnn <= 10000, ]
SSC_parent_US_SV_SNVs$UID <- paste(SSC_parent_US_SV_SNVs$X.Sample, SSC_parent_US_SV_SNVs$X.id, sep='.')
SSC_parent_US_SV_SNVs <- SSC_parent_US_SV_SNVs[!(SSC_parent_US_SV_SNVs$LoF & SSC_parent_US_SV_SNVs$freq_max > 0.01), ]
SSC_parent_US_SV_ExonicSizes <- data.frame(SSC_parent_US_SV_SNVsExonicSizes[[2]])
SSC_parent_US_all_ExonicSizes <- merge(SSC_parent_US_ExonicSizes, SSC_parent_US_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T)
#Combine
SSC_parent_US_all_ExonicSizes <- merge(SSC_parent_US_ExonicSizes, SSC_parent_US_SV_ExonicSizes, by.x='Sample.ID',by.y='Sample.ID',all.x=T, all.y=T)
SSC_parent_US_all_ExonicSizes$exonicSize.x[is.na(SSC_parent_US_all_ExonicSizes$exonicSize.x)] <- 0
SSC_parent_US_all_ExonicSizes$exonicSize.y[is.na(SSC_parent_US_all_ExonicSizes$exonicSize.y)] <- 0
SSC_parent_US_all_ExonicSizes$exonicSize <- pmax(SSC_parent_US_all_ExonicSizes$exonicSize.x, SSC_parent_US_all_ExonicSizes$exonicSize.y)
SSC_parent_US_combined_SNVs <- rbind(SSC_parent_US_SNVs, SSC_parent_US_SV_SNVs)
SSC_parent_US_proc_SNVs <- dplyr::distinct(SSC_parent_US_combined_SNVs, UID, .keep_all=T)


# Filter out FP SNV calls from SV (homozygous dels or < 50bp distance SNV)
SSC_parent_US_SVs <- data.frame(SSC_parent_US_SV_SNVsExonicSizes[[3]])
SSC_parent_US_SVs <- SSC_parent_US_SVs[SSC_parent_US_SVs$length <= 10000
                                       & SSC_parent_US_SVs$Sample.ID %in% SSC_US_meta$Sample.ID, ]

SSC_US_homozyg_SVs <- SSC_parent_US_SVs[SSC_parent_US_SVs$MANTA_GT=='1/1:1',c('sample','CHROM','START','END')]
SSC_US_homozyg_SVs$SVUID <- with(SSC_US_homozyg_SVs, paste0(sample,CHROM,START,END))

SSC_parent_US_proc_SNVs <- SSC_parent_US_proc_SNVs[!with(SSC_parent_US_proc_SNVs,paste0(X.Sample,CHROM,cnvSTARTAnn,cnvENDAnn)) %in% SSC_US_homozyg_SVs$SVUID,]

SSC_parent_US_proc_SNVs$Min_Dist <- by(SSC_parent_US_proc_SNVs, seq_len(nrow(SSC_parent_US_proc_SNVs)),
                                            function(r) r$MIN = min(abs(r$POS - SSC_parent_US_proc_SNVs$POS[SSC_parent_US_proc_SNVs$X.Sample == r$X.Sample & 
                                                                                                              SSC_parent_US_proc_SNVs$CHROM == r$CHROM &
                                                                                                              SSC_parent_US_proc_SNVs$POS != r$POS])) > 50)
SSC_parent_US_proc_SNVs <- SSC_parent_US_proc_SNVs[SSC_parent_US_proc_SNVs$Min_Dist,]

SSC_parent_US_clogit_res <- Get_Log_Reg(SSC_US_meta,
                                        SSC_parent_US_proc_SNVs,
                                        gene_set_data_path,
                                        SSC_parent_US_all_ExonicSizes,
                                        ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                        mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                        type='clogit')
### SSC Proband-Unaffected Siblings
SSC_sibling_metadata_path <- "SSC_metadata.tsv"
SSC_sibling_meta <- Get_Filtered_Metadata(SSC_sibling_metadata_path, child_relation = 'proband',
                                          comparison = 'child-child')
SSC_proband_US_proc_SNVs <- rbind(SSC_parent_US_proc_SNVs[SSC_parent_US_proc_SNVs$X.Sample %in% SSC_sibling_meta$Sample.ID, ],
                                  SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$X.Sample %in% SSC_sibling_meta$Sample.ID, ])
SSC_proband_US_all_ExonicSizes <- rbind(SSC_parent_US_all_ExonicSizes[SSC_parent_US_all_ExonicSizes$Sample.ID %in% SSC_sibling_meta$Sample.ID, ],
                                        SSC_parent_proband_all_ExonicSizes[SSC_parent_proband_all_ExonicSizes$Sample.ID %in% SSC_sibling_meta$Sample.ID, ])

SSC_proband_US_clogit_res <- Get_Log_Reg(SSC_sibling_meta,
                                         SSC_proband_US_proc_SNVs,
                                         gene_set_data_path,
                                         SSC_proband_US_all_ExonicSizes,
                                         snv_count_groups = 'proband-unaffected sibling',
                                         ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                         mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                         type='clogit')

#### CH Event Limited to 0.1% *NOTE: Overwrites previous vars*####
### MSSNG Parent-Proband
MSSNG_parent_proband_0.1CH_SNVs <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$CHfreq <= 0.001, ]
MSSNG_parent_proband_0.1CH_clogit_res <- Get_Log_Reg(MSSNG_meta,
                                                     MSSNG_parent_proband_0.1CH_SNVs,
                                                     gene_set_data_path,
                                                     MSSNG_parent_proband_all_ExonicSizes,
                                                     ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                     mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                     type='clogit')

MSSNG_parent_proband_0.01CH_SNVs <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$CHfreq <= 0.0001, ]
MSSNG_parent_proband_0.01CH_clogit_res <- Get_Log_Reg(MSSNG_meta,
                                                      MSSNG_parent_proband_0.01CH_SNVs,
                                                      gene_set_data_path,
                                                      MSSNG_parent_proband_all_ExonicSizes,
                                                      ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                      mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                      type='clogit')
#SPX test
MSSNG_parent_proband_0.1CH_SNVs_SPX <- MSSNG_parent_proband_proc_SNVs_SPX[MSSNG_parent_proband_proc_SNVs_SPX$CHfreq <= 0.001, ]
MSSNG_parent_proband_0.1CH_clogit_res_SPX <- Get_Log_Reg(MSSNG_meta_SPX,
                                                         MSSNG_parent_proband_0.1CH_SNVs_SPX,
                                                         gene_set_data_path,
                                                         MSSNG_parent_proband_all_ExonicSizes,
                                                         ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                         mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                         type='clogit')

MSSNG_parent_proband_0.01CH_SNVs_SPX <- MSSNG_parent_proband_proc_SNVs_SPX[MSSNG_parent_proband_proc_SNVs_SPX$CHfreq <= 0.0001, ]
MSSNG_parent_proband_0.01CH_clogit_res_SPX <- Get_Log_Reg(MSSNG_meta_SPX,
                                                          MSSNG_parent_proband_0.01CH_SNVs_SPX,
                                                          gene_set_data_path,
                                                          MSSNG_parent_proband_all_ExonicSizes,
                                                          ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                          mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                          type='clogit')

MSSNG_parent_proband_0.1CH_SNVs_MPX <- MSSNG_parent_proband_proc_SNVs_MPX[MSSNG_parent_proband_proc_SNVs_MPX$CHfreq <= 0.001, ]
MSSNG_parent_proband_0.1CH_clogit_res_MPX <- Get_Log_Reg(MSSNG_meta_MPX,
                                                         MSSNG_parent_proband_0.1CH_SNVs_MPX,
                                                         gene_set_data_path,
                                                         MSSNG_parent_proband_all_ExonicSizes,
                                                         ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                         mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                         type='clogit')

MSSNG_parent_proband_0.01CH_SNVs_MPX <- MSSNG_parent_proband_proc_SNVs_MPX[MSSNG_parent_proband_proc_SNVs_MPX$CHfreq <= 0.0001, ]
MSSNG_parent_proband_0.01CH_clogit_res_MPX <- Get_Log_Reg(MSSNG_meta_MPX,
                                                          MSSNG_parent_proband_0.01CH_SNVs_MPX,
                                                          gene_set_data_path,
                                                          MSSNG_parent_proband_all_ExonicSizes,
                                                          ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                          mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                          type='clogit')
### SSC Parent-Proband
SSC_parent_proband_0.1CH_SNVs <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$CHfreq <= 0.001, ]
SSC_parent_proband_0.1CH_clogit_res <- Get_Log_Reg(SSC_meta,
                                                   SSC_parent_proband_0.1CH_SNVs,
                                                   gene_set_data_path,
                                                   SSC_parent_proband_all_ExonicSizes,
                                                   ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                   mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                   type='clogit')

SSC_parent_proband_0.01CH_SNVs <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$CHfreq <= 0.0001, ]
SSC_parent_proband_0.01CH_clogit_res <- Get_Log_Reg(SSC_meta,
                                                 SSC_parent_proband_0.01CH_SNVs,
                                                 gene_set_data_path,
                                                 SSC_parent_proband_all_ExonicSizes,
                                                 ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                 mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                 type='clogit')
### SSC Parent-Unaffected Siblings
SSC_parent_US_0.1CH_SNVs <- SSC_parent_US_proc_SNVs[SSC_parent_US_proc_SNVs$CHfreq <= 0.001, ]
SSC_parent_US_0.1CH_clogit_res <- Get_Log_Reg(SSC_US_meta,
                                              SSC_parent_US_0.1CH_SNVs,
                                              gene_set_data_path,
                                              SSC_parent_US_all_ExonicSizes,
                                              ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                              mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                              type='clogit')

SSC_parent_US_0.01CH_SNVs <- SSC_parent_US_proc_SNVs[SSC_parent_US_proc_SNVs$CHfreq <= 0.0001, ]
SSC_parent_US_0.01CH_clogit_res <- Get_Log_Reg(SSC_US_meta,
                                               SSC_parent_US_0.01CH_SNVs,
                                               gene_set_data_path,
                                               SSC_parent_US_all_ExonicSizes,
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

### SSC Proband-Unaffected Siblings
SSC_proband_US_0.1CH_SNVs <- SSC_proband_US_proc_SNVs[SSC_proband_US_proc_SNVs$CHfreq <= 0.001, ]
SSC_proband_US_0.1CH_clogit_res <- Get_Log_Reg(SSC_sibling_meta,
                                               SSC_proband_US_0.1CH_SNVs,
                                               gene_set_data_path,
                                               SSC_proband_US_all_ExonicSizes,
                                               snv_count_groups = 'proband-unaffected sibling',
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

SSC_proband_US_0.01CH_SNVs <- SSC_proband_US_proc_SNVs[SSC_proband_US_proc_SNVs$CHfreq <= 0.0001, ]
SSC_proband_US_0.01CH_clogit_res <- Get_Log_Reg(SSC_sibling_meta,
                                                SSC_proband_US_0.01CH_SNVs,
                                                gene_set_data_path,
                                                SSC_proband_US_all_ExonicSizes,
                                                snv_count_groups = 'proband-unaffected sibling',
                                                ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                type='clogit')

#### Different SNV MAF Tests ####
MSSNG_parent_proband_proc_SNVs_10p <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$freq_max <= 0.1,]
MSSNG_parent_proband_SNV10p_clogit_res <- Get_Log_Reg(MSSNG_meta,
                                                      MSSNG_parent_proband_proc_SNVs_10p,
                                                      gene_set_data_path,
                                                      MSSNG_parent_proband_all_ExonicSizes,
                                                      ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                      mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                      type='clogit')
MSSNG_parent_proband_proc_SNVs_1p <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$freq_max <= 0.01,]
MSSNG_parent_proband_SNV1p_clogit_res <- Get_Log_Reg(MSSNG_meta,
                                                     MSSNG_parent_proband_proc_SNVs_1p,
                                                     gene_set_data_path,
                                                     MSSNG_parent_proband_all_ExonicSizes,
                                                     ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                     mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                     type='clogit')

SSC_parent_proband_proc_SNVs_10p <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$freq_max <= 0.1,]
SSC_parent_proband_SNV10p_clogit_res <- Get_Log_Reg(SSC_meta,
                                             SSC_parent_proband_proc_SNVs_10p,
                                             gene_set_data_path,
                                             SSC_parent_proband_all_ExonicSizes,
                                             ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                             mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                             type='clogit')
SSC_parent_proband_proc_SNVs_1p <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$freq_max <= 0.01,]
SSC_parent_proband_SNV1p_clogit_res <- Get_Log_Reg(SSC_meta,
                                                    SSC_parent_proband_proc_SNVs_1p,
                                                    gene_set_data_path,
                                                    SSC_parent_proband_all_ExonicSizes,
                                                    ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                                    mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                                    type='clogit')

SSC_parent_US_proc_SNVs_10p <- SSC_parent_US_proc_SNVs[SSC_parent_US_proc_SNVs$freq_max <= 0.1,]
SSC_parent_US_SNV10p_clogit_res <- Get_Log_Reg(SSC_US_meta,
                                        SSC_parent_US_proc_SNVs_10p,
                                        gene_set_data_path,
                                        SSC_parent_US_all_ExonicSizes,
                                        ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                        mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                        type='clogit')
SSC_parent_US_proc_SNVs_1p <- SSC_parent_US_proc_SNVs[SSC_parent_US_proc_SNVs$freq_max <= 0.01,]
SSC_parent_US_SNV1p_clogit_res <- Get_Log_Reg(SSC_US_meta,
                                               SSC_parent_US_proc_SNVs_1p,
                                               gene_set_data_path,
                                               SSC_parent_US_all_ExonicSizes,
                                               ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                               mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                               type='clogit')

SSC_proband_US_proc_SNVs_10p <- SSC_proband_US_proc_SNVs[SSC_proband_US_proc_SNVs$freq_max <= 0.1,]
SSC_proband_US_SNV10p_clogit_res <- Get_Log_Reg(SSC_sibling_meta,
                                         SSC_proband_US_proc_SNVs_10p,
                                         gene_set_data_path,
                                         SSC_proband_US_all_ExonicSizes,
                                         snv_count_groups = 'proband-unaffected sibling',
                                         ref=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn,
                                         mod=Status ~ strata(Family.ID) + Sex + totalCHEvents_syn + Freq,
                                         type='clogit')
SSC_proband_US_proc_SNVs_1p <- SSC_proband_US_proc_SNVs[SSC_proband_US_proc_SNVs$freq_max <= 0.01,]
SSC_proband_US_SNV1p_clogit_res <- Get_Log_Reg(SSC_sibling_meta,
                                                SSC_proband_US_proc_SNVs_1p,
                                                gene_set_data_path,
                                                SSC_proband_US_all_ExonicSizes,
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

reg_resuts <- list(MSSNG_parent_proband_clogit_res=MSSNG_parent_proband_clogit_res,
                   SSC_parent_proband_clogit_res=SSC_parent_proband_clogit_res,
                   SSC_parent_US_clogit_res=SSC_parent_US_clogit_res,
                   SSC_proband_US_clogit_res=SSC_proband_US_clogit_res,
                   MSSNG_parent_proband_0.1CH_clogit_res=MSSNG_parent_proband_0.1CH_clogit_res,
                   SSC_parent_proband_0.1CH_clogit_res=SSC_parent_proband_0.1CH_clogit_res,
                   SSC_parent_US_0.1CH_clogit_res=SSC_parent_US_0.1CH_clogit_res,
                   SSC_proband_US_0.1CH_clogit_res=SSC_proband_US_0.1CH_clogit_res,
                   MSSNG_parent_proband_0.01CH_clogit_res=MSSNG_parent_proband_0.01CH_clogit_res,
                   SSC_parent_proband_0.01CH_clogit_res=SSC_parent_proband_0.01CH_clogit_res,
                   SSC_parent_US_0.01CH_clogit_res=SSC_parent_US_0.01CH_clogit_res,
                   SSC_proband_US_0.01CH_clogit_res=SSC_proband_US_0.01CH_clogit_res,
                   MSSNG_parent_proband_0.1CH_clogit_res_SPX=MSSNG_parent_proband_0.1CH_clogit_res_SPX,
                   MSSNG_parent_proband_0.01CH_clogit_res_SPX=MSSNG_parent_proband_0.01CH_clogit_res_SPX)

slide_data <- format(Get_Regression_Slide_Data(reg_resuts, include_base_reg=F), 
                     scientific=F)

Get_Regression_Slide_Data(list(MSSNG_parent_proband_0.1CH_clogit_res_MPX=MSSNG_parent_proband_0.1CH_clogit_res_MPX,
                               MSSNG_parent_proband_0.01CH_clogit_res_MPX=MSSNG_parent_proband_0.01CH_clogit_res_MPX))
Get_Regression_Slide_Data(list(MSSNG_parent_proband_SNV10p_clogit_res=MSSNG_parent_proband_SNV10p_clogit_res,
                               MSSNG_parent_proband_SNV1p_clogit_res=MSSNG_parent_proband_SNV1p_clogit_res,
                               SSC_parent_proband_SNV10p_clogit_res=SSC_parent_proband_SNV10p_clogit_res,
                               SSC_parent_proband_SNV1p_clogit_res=SSC_parent_proband_SNV1p_clogit_res,
                               SSC_parent_US_SNV10p_clogit_res=SSC_parent_US_SNV10p_clogit_res,
                               SSC_parent_US_SNV1p_clogit_res=SSC_parent_US_SNV1p_clogit_res,
                               SSC_proband_US_SNV10p_clogit_res=SSC_proband_US_SNV10p_clogit_res,
                               SSC_proband_US_SNV1p_clogit_res=SSC_proband_US_SNV1p_clogit_res))
#SNV Counts for all SV:
SNV_in_SV_only_MSSNG <- MSSNG_parent_proband_proc_SNVs[!MSSNG_parent_proband_proc_SNVs$UID %in% MSSNG_parent_proband_SNVs$UID, ]
all_SV_count_MSSNG <- length(SNV_in_SV_only_MSSNG$X.Sample)
syn_SV_count_MSSNG <- length(SNV_in_SV_only_MSSNG$X.Sample[SNV_in_SV_only_MSSNG$effect_priority == 'synonymous SNV'])
mis_SV_count_MSSNG <- length(SNV_in_SV_only_MSSNG$X.Sample[SNV_in_SV_only_MSSNG$effect_priority == 'nonsynonymous SNV'])
lof_SV_count_MSSNG <- length(SNV_in_SV_only_MSSNG$X.Sample[which(SNV_in_SV_only_MSSNG$LoF)])

ssc_all_SNVs_combined <- rbind(SSC_parent_proband_proc_SNVs, SSC_parent_US_proc_SNVs)
ssc_all_SNVs_combined_unique <- dplyr::distinct(ssc_all_SNVs_combined, UID, .keep_all = T)
ssc_all_ExonicSizes_combined <- rbind(SSC_parent_proband_all_ExonicSizes, SSC_parent_US_all_ExonicSizes)
ssc_all_ExonicSizes_combined_unique <- dplyr::distinct(ssc_all_ExonicSizes_combined, Sample.ID, .keep_all = T)
SNV_in_SV_only_SSC <- ssc_all_SNVs_combined_unique[!ssc_all_SNVs_combined_unique$UID %in% SSC_parent_proband_SNVs$UID &
                                                     !ssc_all_SNVs_combined_unique$UID %in% SSC_parent_US_SNVs$UID, ]
all_SV_count_SSC <- length(SNV_in_SV_only_SSC$X.Sample)
syn_SV_count_SSC <- length(SNV_in_SV_only_SSC$X.Sample[SNV_in_SV_only_SSC$effect_priority == 'synonymous SNV'])
mis_SV_count_SSC <- length(SNV_in_SV_only_SSC$X.Sample[SNV_in_SV_only_SSC$effect_priority == 'nonsynonymous SNV'])
lof_SV_count_SSC <- length(SNV_in_SV_only_SSC$X.Sample[which(SNV_in_SV_only_SSC$LoF)])

#SNV Counts for all CNV
SNV_in_CNV_only_MSSNG <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$UID %in% MSSNG_parent_proband_SNVs$UID, ]
all_CNV_count_MSSNG <- length(SNV_in_CNV_only_MSSNG$X.Sample)
syn_CNV_count_MSSNG <- length(SNV_in_CNV_only_MSSNG$X.Sample[SNV_in_CNV_only_MSSNG$effect_priority == 'synonymous SNV'])
mis_CNV_count_MSSNG <- length(SNV_in_CNV_only_MSSNG$X.Sample[SNV_in_CNV_only_MSSNG$effect_priority == 'nonsynonymous SNV'])
lof_CNV_count_MSSNG <- length(SNV_in_CNV_only_MSSNG$X.Sample[which(SNV_in_CNV_only_MSSNG$LoF)])

SNV_in_CNV_only_SSC <- ssc_all_SNVs_combined_unique[ssc_all_SNVs_combined_unique$UID %in% SSC_parent_proband_SNVs$UID |
                                                     ssc_all_SNVs_combined_unique$UID %in% SSC_parent_US_SNVs$UID, ]
all_CNV_count_SSC <- length(SNV_in_CNV_only_SSC$X.Sample)
syn_CNV_count_SSC <- length(SNV_in_CNV_only_SSC$X.Sample[SNV_in_CNV_only_SSC$effect_priority == 'synonymous SNV'])
mis_CNV_count_SSC <- length(SNV_in_CNV_only_SSC$X.Sample[SNV_in_CNV_only_SSC$effect_priority == 'nonsynonymous SNV'])
lof_CNV_count_SSC <- length(SNV_in_CNV_only_SSC$X.Sample[which(SNV_in_CNV_only_SSC$LoF)])

# Other SNV in CNV/SV count metrics
avg_bp_sv_MSSNG <- summary(MSSNG_parent_proband_all_ExonicSizes$exonicSize[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_MSSNG$X.Sample])
avg_bp_cnv_MSSNG <- summary(MSSNG_parent_proband_all_ExonicSizes$exonicSize[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_MSSNG$X.Sample])
avg_bp_sv_SSC <- summary(ssc_all_ExonicSizes_combined_unique$exonicSize[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample])
avg_bp_cnv_SSC <- summary(ssc_all_ExonicSizes_combined_unique$exonicSize[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample])

snv_per_bp_CNV_MSSNG <- all_CNV_count_MSSNG / sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_CNV_only_MSSNG$X.Sample])
snv_per_bp_SV_MSSNG <- all_SV_count_MSSNG / sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% SNV_in_SV_only_MSSNG$X.Sample])

snv_per_bp_CNV_SSC <- all_CNV_count_SSC/ sum(ssc_all_ExonicSizes_combined_unique$exonicSize[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_CNV_only_SSC$X.Sample])
snv_per_bp_SV_SSC <- all_SV_count_SSC / sum(ssc_all_ExonicSizes_combined_unique$exonicSize[ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SNV_in_SV_only_SSC$X.Sample])

# No. of CNV & SV counts
MSSNG_parent_proband_CNVs_ILMN <- data.frame(MSSNG_parent_proband_SNVsExonicSizes[[3]])
MSSNG_parent_proband_CNVs_ILMN <- MSSNG_parent_proband_CNVs_ILMN[MSSNG_parent_proband_CNVs_ILMN$freq_max/100 <= 0.01
                                                                 & MSSNG_parent_proband_CNVs_ILMN$Sample.ID %in% MSSNG_meta$Sample.ID, ]
MSSNG_parent_proband_CNVs_CG <- data.frame(MSSNG_parent_proband_SNVsExonicSizes[[4]])
MSSNG_parent_proband_CNVs_CG <- MSSNG_parent_proband_CNVs_CG[MSSNG_parent_proband_CNVs_CG$freq_max/100 <= 0.01
                                                             & MSSNG_parent_proband_CNVs_CG$Sample.ID %in% MSSNG_meta$Sample.ID, ]
SSC_parent_proband_CNVs <- data.frame(SSC_parent_proband_SNVsExonicSizes[[3]])
SSC_parent_proband_CNVs <- SSC_parent_proband_CNVs[SSC_parent_proband_CNVs$freq_max/100 <= 0.01
                                                   & SSC_parent_proband_CNVs$Sample.ID %in% SSC_meta$Sample.ID, ]
SSC_parent_US_CNVs <- data.frame(SSC_parent_US_SNVsExonicSizes[[3]])
SSC_parent_US_CNVs <- SSC_parent_US_CNVs[SSC_parent_US_CNVs$freq_max/100 <= 0.01
                                         & SSC_parent_US_CNVs$Sample.ID %in% SSC_US_meta$Sample.ID, ]
SSC_parent_child_CNVs <- rbind(SSC_parent_proband_CNVs, SSC_parent_US_CNVs)
SSC_parent_child_CNVs <- SSC_parent_child_CNVs[!duplicated(paste(SSC_parent_child_CNVs$Sample.ID, 
                                                            SSC_parent_child_CNVs$CHROM,
                                                            SSC_parent_child_CNVs$START,
                                                            SSC_parent_child_CNVs$END, sep='.')), ]


MSSNG_father_IDs <- MSSNG_meta$Sample.ID[MSSNG_meta$Relation == 'father']
MSSNG_mother_IDs <- MSSNG_meta$Sample.ID[MSSNG_meta$Relation == 'mother']
MSSNG_proband_IDs <- MSSNG_meta$Sample.ID[MSSNG_meta$Relation == 'proband']
MSSNG_CNV_count <- nrow(MSSNG_parent_proband_CNVs_ILMN) + nrow(MSSNG_parent_proband_CNVs_CG)
MSSNG_CNV_father_ILMN <- MSSNG_parent_proband_CNVs_ILMN[MSSNG_parent_proband_CNVs_ILMN$Sample.ID %in% MSSNG_father_IDs, ]
MSSNG_CNV_father_CG <- MSSNG_parent_proband_CNVs_CG[MSSNG_parent_proband_CNVs_CG$Sample.ID %in% MSSNG_father_IDs, ]
MSSNG_CNV_count_father <- nrow(MSSNG_CNV_father_ILMN) + nrow(MSSNG_CNV_father_CG)
MSSNG_CNV_mother_ILMN <- MSSNG_parent_proband_CNVs_ILMN[MSSNG_parent_proband_CNVs_ILMN$Sample.ID %in% MSSNG_mother_IDs, ]
MSSNG_CNV_mother_CG <- MSSNG_parent_proband_CNVs_CG[MSSNG_parent_proband_CNVs_CG$Sample.ID %in% MSSNG_mother_IDs, ]
MSSNG_CNV_count_mother <- nrow(MSSNG_CNV_mother_ILMN) + nrow(MSSNG_CNV_mother_CG)
MSSNG_CNV_proband_ILMN <- MSSNG_parent_proband_CNVs_ILMN[MSSNG_parent_proband_CNVs_ILMN$Sample.ID %in% MSSNG_proband_IDs, ]
MSSNG_CNV_proband_CG <- MSSNG_parent_proband_CNVs_CG[MSSNG_parent_proband_CNVs_CG$Sample.ID %in% MSSNG_proband_IDs, ]
MSSNG_CNV_count_proband <- nrow(MSSNG_CNV_proband_ILMN) + nrow(MSSNG_CNV_proband_CG)

MSSNG_SV_count <- nrow(MSSNG_parent_proband_SVs_ILMN) + nrow(MSSNG_parent_proband_SVs_CG)
MSSNG_SV_father_ILMN <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$Sample.ID %in% MSSNG_father_IDs, ]
MSSNG_SV_father_CG <- MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$Sample.ID %in% MSSNG_father_IDs, ]
MSSNG_SV_count_father <- nrow(MSSNG_SV_father_ILMN) + nrow(MSSNG_SV_father_CG)
MSSNG_SV_mother_ILMN <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$Sample.ID %in% MSSNG_mother_IDs, ]
MSSNG_SV_mother_CG <- MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$Sample.ID %in% MSSNG_mother_IDs, ]
MSSNG_SV_count_mother <- nrow(MSSNG_SV_mother_ILMN) + nrow(MSSNG_SV_mother_CG)
MSSNG_SV_proband_ILMN <- MSSNG_parent_proband_SVs_ILMN[MSSNG_parent_proband_SVs_ILMN$Sample.ID %in% MSSNG_proband_IDs, ]
MSSNG_SV_proband_CG <- MSSNG_parent_proband_SVs_CG[MSSNG_parent_proband_SVs_CG$Sample.ID %in% MSSNG_proband_IDs, ]
MSSNG_SV_count_proband <- nrow(MSSNG_SV_proband_ILMN) + nrow(MSSNG_SV_proband_CG)

SSC_parent_child_SVs <- rbind(SSC_parent_proband_SVs, SSC_parent_US_SVs)
SSC_parent_child_SVs <- SSC_parent_child_SVs[!duplicated(paste(SSC_parent_child_SVs$Sample.ID, 
                                                               SSC_parent_child_SVs$CHROM,
                                                               SSC_parent_child_SVs$START,
                                                               SSC_parent_child_SVs$END, sep='.')), ]



SSC_father_IDs <- SSC_meta$Sample.ID[SSC_meta$Relation == 'father']
SSC_mother_IDs <- SSC_meta$Sample.ID[SSC_meta$Relation == 'mother']
SSC_proband_IDs <- SSC_meta$Sample.ID[SSC_meta$Relation == 'proband']
SSC_father_US_IDs <- SSC_US_meta$Sample.ID[SSC_US_meta$Relation == 'father']
SSC_mother_US_IDs <- SSC_US_meta$Sample.ID[SSC_US_meta$Relation == 'mother']
SSC_US_IDs <- SSC_US_meta$Sample.ID[!SSC_US_meta$Relation %in% c('father', 'mother')]

SSC_CNV_count <- nrow(SSC_parent_child_CNVs)
SSC_CNV_father <- SSC_parent_proband_CNVs[SSC_parent_proband_CNVs$Sample.ID %in% SSC_father_IDs, ]
SSC_CNV_count_father <- nrow(SSC_CNV_father)
SSC_CNV_mother <- SSC_parent_proband_CNVs[SSC_parent_proband_CNVs$Sample.ID %in% SSC_mother_IDs, ]
SSC_CNV_count_mother <- nrow(SSC_CNV_mother)
SSC_CNV_proband <- SSC_parent_proband_CNVs[SSC_parent_proband_CNVs$Sample.ID %in% SSC_proband_IDs, ]
SSC_CNV_count_proband <- nrow(SSC_CNV_proband)
SSC_CNV_father_US <- SSC_parent_US_CNVs[SSC_parent_US_CNVs$Sample.ID %in% SSC_father_US_IDs, ]
SSC_CNV_count_father_US <- nrow(SSC_CNV_father_US)
SSC_CNV_mother_US <- SSC_parent_US_CNVs[SSC_parent_US_CNVs$Sample.ID %in% SSC_mother_US_IDs, ]
SSC_CNV_count_mother_US <- nrow(SSC_CNV_mother_US)
SSC_CNV_US <- SSC_parent_US_CNVs[SSC_parent_US_CNVs$Sample.ID %in% SSC_US_IDs, ]
SSC_CNV_count_US <- nrow(SSC_CNV_US)

SSC_SV_count <- nrow(SSC_parent_child_SVs)
SSC_SV_father <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$Sample.ID %in% SSC_father_IDs, ]
SSC_SV_count_father <- nrow(SSC_SV_father)
SSC_SV_mother <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$Sample.ID %in% SSC_mother_IDs, ]
SSC_SV_count_mother <- nrow(SSC_SV_mother)
SSC_SV_proband <- SSC_parent_proband_SVs[SSC_parent_proband_SVs$Sample.ID %in% SSC_proband_IDs, ]
SSC_SV_count_proband <- nrow(SSC_SV_proband)
SSC_SV_father_US <- SSC_parent_US_SVs[SSC_parent_US_SVs$Sample.ID %in% SSC_father_US_IDs, ]
SSC_SV_count_father_US <- nrow(SSC_SV_father_US)
SSC_SV_mother_US <- SSC_parent_US_SVs[SSC_parent_US_SVs$Sample.ID %in% SSC_mother_US_IDs, ]
SSC_SV_count_mother_US <- nrow(SSC_SV_mother_US)
SSC_SV_US <- SSC_parent_US_SVs[SSC_parent_US_SVs$Sample.ID %in% SSC_US_IDs, ]
SSC_SV_count_US <- nrow(SSC_SV_US)

## Total exon lengths from each individual group
MSSNG_CNV_father_totalExonLen <- sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
  MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_CNV_father_ILMN$Sample.ID]) +
  sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
    MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_CNV_father_CG$Sample.ID])
MSSNG_CNV_mother_totalExonLen <- sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
  MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_CNV_mother_ILMN$Sample.ID]) +
  sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
    MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_CNV_mother_CG$Sample.ID])
MSSNG_CNV_proband_totalExonLen <- sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
  MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_CNV_proband_ILMN$Sample.ID]) +
  sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
    MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_CNV_proband_CG$Sample.ID])

MSSNG_SV_father_totalExonLen <- sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
  MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_SV_father_ILMN$Sample.ID]) +
  sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
    MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_SV_father_CG$Sample.ID])
MSSNG_SV_mother_totalExonLen <- sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
  MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_SV_mother_ILMN$Sample.ID]) +
  sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
    MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_SV_mother_CG$Sample.ID])
MSSNG_SV_proband_totalExonLen <- sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
  MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_SV_proband_ILMN$Sample.ID]) +
  sum(MSSNG_parent_proband_all_ExonicSizes$exonicSize[
    MSSNG_parent_proband_all_ExonicSizes$Sample.ID %in% MSSNG_SV_proband_CG$Sample.ID])


SSC_CNV_father_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_CNV_father$Sample.ID])
SSC_CNV_mother_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_CNV_mother$Sample.ID])
SSC_CNV_proband_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_CNV_proband$Sample.ID])
SSC_CNV_father_US_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_CNV_father_US$Sample.ID])
SSC_CNV_mother_US_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_CNV_mother_US$Sample.ID])
SSC_CNV_US_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_CNV_US$Sample.ID])

SSC_SV_father_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_SV_father$Sample.ID])
SSC_SV_mother_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_SV_mother$Sample.ID])
SSC_SV_proband_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_SV_proband$Sample.ID])
SSC_SV_father_US_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_SV_father_US$Sample.ID])
SSC_SV_mother_US_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_SV_mother_US$Sample.ID])
SSC_SV_US_totalExonLen <- sum(ssc_all_ExonicSizes_combined_unique$exonicSize[
  ssc_all_ExonicSizes_combined_unique$Sample.ID %in% SSC_SV_US$Sample.ID])

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

mssng_parent_proband_indiv_data <- Get_Number_Of_Individual_Comparison(MSSNG_parent_proband_proc_SNVs, 
                                                                       MSSNG_meta,
                                                                       comparison='parent-child', 
                                                                       child_relation='proband')
ssc_parent_proband_indiv_data <- Get_Number_Of_Individual_Comparison(SSC_parent_proband_proc_SNVs, 
                                                                       SSC_meta,
                                                                       comparison='parent-child', 
                                                                       child_relation='proband')
ssc_parent_US_indiv_data <- Get_Number_Of_Individual_Comparison(SSC_parent_US_proc_SNVs, 
                                                                     SSC_US_meta,
                                                                     comparison='parent-child', 
                                                                     child_relation='unaffected sibling')
ssc_proband_US_indiv_data <- Get_Number_Of_Individual_Comparison(SSC_proband_US_proc_SNVs, 
                                                                     SSC_sibling_meta,
                                                                     comparison='child-child', 
                                                                     child_relation='proband')

mssng_parent_proband_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(MSSNG_parent_proband_0.1CH_SNVs,
                                                                       MSSNG_meta,
                                                                       comparison='parent-child', 
                                                                       child_relation='proband')
ssc_parent_proband_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(SSC_parent_proband_0.1CH_SNVs, 
                                                                     SSC_meta,
                                                                     comparison='parent-child', 
                                                                     child_relation='proband')
ssc_parent_US_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(SSC_parent_US_0.1CH_SNVs, 
                                                                SSC_US_meta,
                                                                comparison='parent-child', 
                                                                child_relation='unaffected sibling')
ssc_proband_US_indiv_0.1_data <- Get_Number_Of_Individual_Comparison(SSC_proband_US_0.1CH_SNVs, 
                                                                 SSC_sibling_meta,
                                                                 comparison='child-child', 
                                                                 child_relation='proband')

Generate_Indiv_Count_Bar_Plots(mssng_parent_proband_indiv_data[[1]], 'mssng_parent_proband_indiv_data')
Generate_Indiv_Count_Bar_Plots(ssc_parent_proband_indiv_data[[1]], 'ssc_parent_proband_indiv_data')
Generate_Indiv_Count_Bar_Plots(ssc_parent_US_indiv_data[[1]], 'ssc_parent_US_indiv_data', comparison_type='parent-child-US')
Generate_Indiv_Count_Bar_Plots(ssc_proband_US_indiv_data[[1]], 'ssc_proband_US_indiv_data', comparison_type='child-child')
Generate_Indiv_Count_Bar_Plots(mssng_parent_proband_indiv_0.1_data[[1]], 'mssng_parent_proband_indiv_0.1_data')
Generate_Indiv_Count_Bar_Plots(ssc_parent_proband_indiv_0.1_data[[1]], 'ssc_parent_proband_indiv_0.1_data')
Generate_Indiv_Count_Bar_Plots(ssc_parent_US_indiv_0.1_data[[1]], 'ssc_parent_US_indiv_0.1_data', comparison_type='parent-child-US')
Generate_Indiv_Count_Bar_Plots(ssc_proband_US_indiv_0.1_data[[1]], 'ssc_proband_US_indiv_0.1_data', comparison_type='child-child')

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
gs_reg_result_list <- list(MSSNG_parent_proband_0.1CH_clogit_res=MSSNG_parent_proband_0.1CH_clogit_res,
                           SSC_parent_proband_0.1CH_clogit_res=SSC_parent_proband_0.1CH_clogit_res,
                           SSC_parent_US_0.1CH_clogit_res=SSC_parent_US_0.1CH_clogit_res,
                           SSC_proband_US_0.1CH_clogit_res=SSC_proband_US_0.1CH_clogit_res)

gs_reg_data <- format(Get_Gene_Set_Regression_Data(gs_reg_result_list), scientific=F)

## Genes with multiple hits/multiple samples
MSSNG_gene_summary <- MSSNG_parent_proband_proc_SNVs %>% filter(X.Sample%in%MSSNG_proband_IDs & effect_priority=='nonsynonymous SNV' & gnomAD_pRec >= 0.9 & gnomAD_oe_lof_upper >= 0.35) %>% 
                      group_by(gene_symbol) %>% summarise(Proband=length(unique(MSSNG_parent_proband_proc_SNVs$X.Sample[MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_proband_IDs & MSSNG_parent_proband_proc_SNVs$gene_symbol %in% gene_symbol & MSSNG_parent_proband_proc_SNVs$effect_priority=='nonsynonymous SNV'])),
                                                          Father=length(unique(MSSNG_parent_proband_proc_SNVs$X.Sample[MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_father_IDs & MSSNG_parent_proband_proc_SNVs$gene_symbol %in% gene_symbol & MSSNG_parent_proband_proc_SNVs$effect_priority=='nonsynonymous SNV'])),
                                                          Mother=length(unique(MSSNG_parent_proband_proc_SNVs$X.Sample[MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_mother_IDs & MSSNG_parent_proband_proc_SNVs$gene_symbol %in% gene_symbol & MSSNG_parent_proband_proc_SNVs$effect_priority=='nonsynonymous SNV']))) %>% arrange(-Proband)
MSSNG_gene_summary <- data.frame(MSSNG_gene_summary[MSSNG_gene_summary$Proband > 1,])

SSC_gene_summary <- ssc_all_SNVs_combined_unique %>% filter(X.Sample%in%SSC_proband_IDs & effect_priority=='nonsynonymous SNV' & gnomAD_pRec >= 0.9 & gnomAD_oe_lof_upper >= 0.35) %>% 
  group_by(gene_symbol) %>% summarise(Proband=length(unique(ssc_all_SNVs_combined_unique$X.Sample[ssc_all_SNVs_combined_unique$X.Sample %in% SSC_proband_IDs & ssc_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & ssc_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV'])), 
                                      Father=length(unique(ssc_all_SNVs_combined_unique$X.Sample[ssc_all_SNVs_combined_unique$X.Sample %in% SSC_father_IDs & ssc_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & ssc_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV'])),
                                      Mother=length(unique(ssc_all_SNVs_combined_unique$X.Sample[ssc_all_SNVs_combined_unique$X.Sample %in% SSC_mother_IDs & ssc_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & ssc_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV'])),
                                      US=length(unique(ssc_all_SNVs_combined_unique$X.Sample[ssc_all_SNVs_combined_unique$X.Sample %in% SSC_US_IDs & ssc_all_SNVs_combined_unique$gene_symbol %in% gene_symbol & ssc_all_SNVs_combined_unique$effect_priority=='nonsynonymous SNV']))) %>% arrange(-Proband)
SSC_gene_summary <- data.frame(SSC_gene_summary[SSC_gene_summary$Proband > 1,])
write.csv(MSSNG_gene_summary, '../DT/MSSNG_proband_genes_missense_summary.csv')
write.csv(SSC_gene_summary, '../DT/SSC_proband_genes_missense_summary.csv')

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
# 
Find_Gene_Count_DT <- function(full_table, IDs, gene) {
  full_table <- full_table[full_table$gnomAD_pRec >= 0.9 & 
                             full_table$gnomAD_oe_lof_upper >= 0.35,]
  gene_count <- length(full_table$CHROM[full_table$X.Sample %in% IDs 
                                        & full_table$gene_symbol %in% gene])
  
  return (gene_count)
}


MSSNG_proband_nonsyn_genes_DT <- MSSNG_parent_proband_proc_SNVs[MSSNG_parent_proband_proc_SNVs$effect_priority == 'nonsynonymous SNV' & 
                                                                  MSSNG_parent_proband_proc_SNVs$X.Sample %in% MSSNG_proband_IDs &
                                                                  MSSNG_parent_proband_proc_SNVs$gnomAD_pRec >= 0.9 & 
                                                                  MSSNG_parent_proband_proc_SNVs$gnomAD_oe_lof_upper >= 0.35,
                                               c('X.Sample','gene_symbol')]
MSSNG_proband_nonsyn_genes_DT <- data.frame(MSSNG_proband_nonsyn_genes_DT %>% 
                                              distinct(.keep_all=T) %>%
                                              group_by(gene_symbol) %>% 
                                              filter(n() > 1) %>%
                                              summarize(Freq=n(), 
                                                        Father=Find_Gene_Count_DT(MSSNG_parent_proband_proc_SNVs,MSSNG_father_IDs,gene_symbol),
                                                        Mother=Find_Gene_Count_DT(MSSNG_parent_proband_proc_SNVs,MSSNG_mother_IDs,gene_symbol)) %>% 
                                              arrange(-Freq))
write.table(MSSNG_proband_nonsyn_genes_DT, "../DT/GeneCounts/MSSNG_proband_nonsyn_genes_DT.tsv", sep="\t", row.names=F, quote=F, col.names=T)

SSC_proband_nonsyn_genes_DT <- SSC_parent_proband_proc_SNVs[SSC_parent_proband_proc_SNVs$effect_priority == 'nonsynonymous SNV' & 
                                                                  SSC_parent_proband_proc_SNVs$X.Sample %in% SSC_proband_IDs &
                                                                  SSC_parent_proband_proc_SNVs$gnomAD_pRec >= 0.9 & 
                                                                  SSC_parent_proband_proc_SNVs$gnomAD_oe_lof_upper >= 0.35,
                                                                c('X.Sample','gene_symbol')]
SSC_proband_nonsyn_genes_DT <- data.frame(SSC_proband_nonsyn_genes_DT %>% 
                                              distinct(.keep_all=T) %>%
                                              group_by(gene_symbol) %>% 
                                              filter(n() > 1) %>%
                                              summarize(Freq=n(), 
                                                        Father=Find_Gene_Count_DT(SSC_parent_proband_proc_SNVs,SSC_father_IDs,gene_symbol),
                                                        Mother=Find_Gene_Count_DT(SSC_parent_proband_proc_SNVs,SSC_mother_IDs,gene_symbol),
                                                        US=Find_Gene_Count_DT(SSC_parent_US_proc_SNVs,SSC_US_IDs,gene_symbol)) %>% 
                                              arrange(-Freq))
write.table(SSC_proband_nonsyn_genes_DT, "../DT/GeneCounts/SSC_proband_nonsyn_genes_DT.tsv", sep="\t", row.names=F, quote=F, col.names=T)


