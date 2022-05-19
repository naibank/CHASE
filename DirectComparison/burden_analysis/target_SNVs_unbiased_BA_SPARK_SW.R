#### SPARK burden analysis extract SNVs with:
# - SNV freq 100%
# - CNV freq < 1%
# - all variants
# - parent or proband SNVs 

library(yaml)
library(survival)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/Users/shaniawu/Desktop/CHASE/SPARK_BA/ExtractSNVs/")

## from Faraz code:
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
  Subjects_to_exclude <- fread('../Data/excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) 
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('../Data/samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) 
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("../Data/MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("../Data/MSSNG+SSC.CNVs.tsv")
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
    # Get only samples from families with at least 1 of each sibling
    proband_fam_IDs <- meta$Family.ID[meta$Relation == 'proband']
    sibling_fam_IDs <- meta$Family.ID[meta$Relation %in% c('unaffected sibling', 'other sibling')]
    meta <- meta[meta$Family.ID %in% proband_fam_IDs &
                   meta$Family.ID %in% sibling_fam_IDs,]
  }
  
  meta$Status <- ifelse(meta$Relation %in% child_relations, 1, 0)
  return (meta)
}


Get_Target_SNVs <- function(meta, snvs, comparison, SNV.freq) {

  ## filter SNVs (from Faraz Get_Log_Reg):
  snvs <- snvs[snvs$X.Sample %in% meta$Sample.ID,]
  #snvs <- snvs[which(snvs$freq_max <= SNV.freq),]
  #snvs <- snvs[!duplicated(snvs[, "UID"]),]
  snvs <- snvs[!duplicated(snvs[, c("X.Sample", "X.id")]),]
  snvs <- snvs[which(snvs$gnomAD_pRec >= 0.9 & snvs$gnomAD_oe_lof_upper >= 0.35), ]

  ## get SNVs
  if (comparison == 'parent-proband') {
    proband_IDs <- meta$Sample.ID[!meta$Relation == 'mother' & !meta$Relation == 'father']
    father_IDs <- meta$Sample.ID[meta$Relation == 'father']
    mother_IDs <- meta$Sample.ID[meta$Relation == 'mother']

    proband_SNVs <- snvs[which(snvs$X.Sample %in% proband_IDs),]
    proband_SNVs$Relation <- "proband"
    father_SNVs <- snvs[which(snvs$X.Sample %in% father_IDs),]
    father_SNVs$Relation <- "father"
    mother_SNVs <- snvs[which(snvs$X.Sample %in% mother_IDs),]
    mother_SNVs$Relation <- "mother"
    
    parent.proband.SNVs <- rbind(proband_SNVs, father_SNVs, mother_SNVs)

    write.csv(parent.proband.SNVs, "SPARK.parent.proband.SNVs.SNV1.CNV0.01.csv", row.names=F)
    # write.table(proband_SNVs, "SPARK.proband.SNVs.SNV1.CNV0.01.tsv", row.names=F)
    # write.table(father_SNVs, "SPARK.father.SNVs.SNV1.CNV0.01.tsv", row.names=F)
    # write.table(mother_SNVs, "SPARK.mother.SNVs.SNV1.CNV0.01.tsv", row.names=F)
  }
  else if (comparison == 'parent-unaffected sibling') {
    father_IDs <- meta$Sample.ID[meta$Relation == 'father']
    mother_IDs <- meta$Sample.ID[meta$Relation == 'mother']
    US_IDs <- meta$Sample.ID[meta$Relation == 'unaffected sibling']

    father_SNVs <- snvs[which(snvs$X.Sample %in% father_IDs),]
    father_SNVs$Relation <- "father"
    mother_SNVs <- snvs[which(snvs$X.Sample %in% mother_IDs),]
    mother_SNVs$Relation <- "mother"
    US_SNVs <-  snvs[which(snvs$X.Sample %in% US_IDs),]
    US_SNVs$Relation <- "unaffected sibling"
    
    parent.US.SNVs <- rbind(father_SNVs, mother_SNVs, US_SNVs)
    write.csv(parent.US.SNVs, "SPARK.parent.US.SNVs.SNV1.CNV0.01.csv", row.names=F)
    
    # write.table(US_SNVs, "SPARK.US.SNVs.SNV1.CNV0.01.tsv", row.names=F)
    # write.table(father_SNVs, "SPARK.father.SNVs.SNV1.CNV0.01.tsv", row.names=F)
    # write.table(mother_SNVs, "SPARK.mother.SNVs.SNV1.CNV0.01.tsv", row.names=F)
  }
}


###### Import meta, SNV, CNV data from SNVsExonicSizes.yaml 

### SPARK Parent-Probands
SPARK_meta <- data.frame()
SPARK_parent_proband_SNVs <- data.frame()
SPARK_parent_proband_ExonicSizes <- data.frame()
SPARK_parent_proband_CNVs <- data.frame()

for (i in 1:3) {
  SPARK_metadata_path <- paste("../Data/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  SPARK_meta <- rbind(SPARK_meta, Get_Filtered_Metadata(SPARK_metadata_path, child_relation = 'proband'))
  
  SPARK_SNVsExonicSizes_path <- paste('../Data/SPARKWGS',i,'_parent_proband_SNVsExonicSizes.yaml',sep='')
  SPARK_parent_proband_SNVsExonicSizes_curr <- yaml::yaml.load_file(SPARK_SNVsExonicSizes_path)
  SPARK_parent_proband_SNVs <- rbind(SPARK_parent_proband_SNVs,
                                     data.frame(SPARK_parent_proband_SNVsExonicSizes_curr[[1]]))
  SPARK_parent_proband_ExonicSizes <- rbind(SPARK_parent_proband_ExonicSizes,
                                            data.frame(SPARK_parent_proband_SNVsExonicSizes_curr[[2]]))
  SPARK_parent_proband_CNVs <- rbind(SPARK_parent_proband_CNVs,
                                     data.frame(SPARK_parent_proband_SNVsExonicSizes_curr[[3]]))
}

SPARK_parent_proband_CNVs$freq_max <- Determine_CNV_MAF(SPARK_parent_proband_CNVs) # add freq_max col to CNVs
SPARK_parent_proband_SNVs <- SPARK_parent_proband_SNVs[
  SPARK_parent_proband_SNVs$effect_priority %in% c('synonymous SNV','nonsynonymous SNV')
                                                       | SPARK_parent_proband_SNVs$LoF,] # get only nonsyn, syn, or LoF SNVs

SPARK_parent_proband_SNVs$UID <- paste(SPARK_parent_proband_SNVs$X.Sample, SPARK_parent_proband_SNVs$X.id, sep='.')  # UID = ID.location-mutation; unique for each SNV
SPARK_parent_proband_SNVs <- SPARK_parent_proband_SNVs[!(SPARK_parent_proband_SNVs$LoF & SPARK_parent_proband_SNVs$freq_max > 0.01), ] # filter out LoF SNVs with freq > 1%
SPARK_parent_proband_SNVs <- SPARK_parent_proband_SNVs[!is.na(SPARK_parent_proband_SNVs$X.Sample),]

## Filter out false positive SNVs calls (those <50bp away from each other):
SPARK_parent_proband_proc_SNVs <- SPARK_parent_proband_SNVs  
SPARK_parent_proband_proc_SNVs$Min_Dist <- by(SPARK_parent_proband_proc_SNVs, seq_len(nrow(SPARK_parent_proband_proc_SNVs)),
                                              function(r) r$MIN = min(abs(r$POS - SPARK_parent_proband_proc_SNVs$POS[
                                                SPARK_parent_proband_proc_SNVs$X.Sample == r$X.Sample &
                                                  SPARK_parent_proband_proc_SNVs$CHROM == r$CHROM &
                                                  SPARK_parent_proband_proc_SNVs$POS != r$POS])) > 50)
SPARK_parent_proband_proc_SNVs <- SPARK_parent_proband_proc_SNVs[SPARK_parent_proband_proc_SNVs$Min_Dist,]
 
## Make proband, father, mother, target SNV tables
Get_Target_SNVs(SPARK_meta, SPARK_parent_proband_proc_SNVs, comparison = "parent-proband", SNV.freq = 1)



### Check SNV counts with plot counts
 SNV.count.table <- data.frame()

 file.path <- "SPARK.parent.proband.SNVs.SNV1.CNV0.01.csv"
 snv.table <- fread(file.path, data.table=F)

for (relation in c("proband", "father", "mother")){
 snv.table.relation <- snv.table[which(snv.table$Relation == relation),]
 message(nrow(snv.table.relation[!duplicated(snv.table.relation$X.Sample)]))
 all.count <- nrow(snv.table.relation)
 nonsyn.count <- nrow(snv.table.relation[which(snv.table.relation$effect_priority == "nonsynonymous SNV"),])
 syn.count <- nrow(snv.table.relation[which(snv.table.relation$effect_priority == "synonymous SNV"),])
 count.table.row <- data.frame("relation" = relation, "all variants" = all.count, "nonsynonymous" = nonsyn.count,
                            "synonymous" = syn.count)
 SNV.count.table <- rbind(SNV.count.table, count.table.row)
}

write.csv(SNV.count.table, "parent.proband.SNV.count.table.csv", row.names=F)

SNV.count.table <- data.frame()

file.path <- "SPARK.parent.US.SNVs.SNV1.CNV0.01.csv"
snv.table <- fread(file.path, data.table=F)

for (relation in c("unaffected sibling", "father", "mother")){
  snv.table.relation <- snv.table[which(snv.table$Relation == relation),]
  message(nrow(snv.table.relation))
  all.count <- nrow(snv.table.relation)
  nonsyn.count <- nrow(snv.table.relation[which(snv.table.relation$effect_priority == "nonsynonymous SNV"),])
  syn.count <- nrow(snv.table.relation[which(snv.table.relation$effect_priority == "synonymous SNV"),])
  count.table.row <- data.frame("relation" = relation, "all variants" = all.count, "nonsynonymous" = nonsyn.count,
                                "synonymous" = syn.count)
  SNV.count.table <- rbind(SNV.count.table, count.table.row)
}

write.csv(SNV.count.table, "parent.US.SNV.count.table.csv", row.names=F)


### SPARK Parent-Unaffected Siblings
SPARK_US_meta <- data.frame()
SPARK_parent_US_SNVs <- data.frame()
SPARK_parent_US_ExonicSizes <- data.frame()
SPARK_parent_US_CNVs <- data.frame()
for (i in 1:3) {
  SPARK_US_metadata_path <- paste("../Data/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  SPARK_US_meta <- rbind(SPARK_US_meta, Get_Filtered_Metadata(SPARK_US_metadata_path, child_relation = 'unaffected sibling'))

  SPARK_US_SNVsExonicSizes_path <- paste('../Data/SPARKWGS',i,'_parent_US_SNVsExonicSizes.yaml',sep='')
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

SPARK_parent_US_proc_SNVs <- SPARK_parent_US_SNVs
SPARK_parent_US_proc_SNVs$Min_Dist <- by(SPARK_parent_US_proc_SNVs, seq_len(nrow(SPARK_parent_US_proc_SNVs)),
                                         function(r) r$MIN = min(abs(r$POS - SPARK_parent_US_proc_SNVs$POS[
                                           SPARK_parent_US_proc_SNVs$X.Sample == r$X.Sample & 
                                             SPARK_parent_US_proc_SNVs$CHROM == r$CHROM &
                                             SPARK_parent_US_proc_SNVs$POS != r$POS])) > 50)
SPARK_parent_US_proc_SNVs <- SPARK_parent_US_proc_SNVs[SPARK_parent_US_proc_SNVs$Min_Dist,]


Get_Target_SNVs(SPARK_US_meta, SPARK_parent_US_proc_SNVs, comparison = "parent-unaffected sibling", SNV.freq = 1)

