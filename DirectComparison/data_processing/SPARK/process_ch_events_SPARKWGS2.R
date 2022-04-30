library(yaml)
library(survival)
library(GenomicRanges)
library(Repitools)
library(data.table)
library(dplyr)

setwd("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis")

Get_Failed_QC_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('./Data/excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('./Data/samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("./Data/MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("./Data/MSSNG+SSC.CNVs.tsv")
  ls_exclude <- unique(rbind(ls1, ls2))
  
  full_metadata <- fread(file, data.table = F)
  familyID_exclude <- full_metadata$`Family ID`[which(full_metadata$`Sample ID` %in% ls_exclude$Sample)]
  sampleID_exclude <- full_metadata$`Sample ID`[which(full_metadata$`Family ID` %in% familyID_exclude)]
  
  return (sampleID_exclude)
}

Filter_ASD_Associated_CNVs_and_LOFs <- function(df, ds='MSSNG') {
  if (ds == "MSSNG") {
    file <- "./Data/MSSNG_metadata.tsv"
  }
  else if (ds == "SSC") {
    file <- "./Data/SSC_metadata.tsv"
  }
  sampleID_exclude <- Get_SampleID_ASD_Associated_CNVs_and_LOFs(file)
  df <- df[which(!df$`X.Sample` %in% sampleID_exclude),]
  
  return (df)
}

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

Determine_SV_MAF <- function(SVs, is_CG=F) {
  if (is_CG) {
    SV_MAF <- pmax(SVs$CGparentalPercFreq_50percRecipOverlap,
                   pmax(SVs$otgMantaPercFreq_90percRecipOverlap,
                        pmax(SVs$otgDellyPercFreq_90percRecipOverlap,
                             pmax(SVs$sscMantaPercFreq_90percRecipOverlap,
                                  pmax(SVs$sscDellyPercFreq_90percRecipOverlap,
                                       pmax(SVs$svMantaXPercFreq_90percRecipOverlap,
                                            pmax(SVs$svManta2PercFreq_90percRecipOverlap,
                                                 pmax(SVs$svDellyXPercFreq_90percRecipOverlap,
                                                      SVs$svDelly2PercFreq_90percRecipOverlap, na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T)/100
  }
  else {
    SV_MAF <- pmax(SVs$CGparentalPercFreq_90percRecipOverlap,
                   pmax(SVs$otgMantaPercFreq_90percRecipOverlap,
                        pmax(SVs$otgDellyPercFreq_90percRecipOverlap,
                             pmax(SVs$sscMantaPercFreq_90percRecipOverlap,
                                  pmax(SVs$sscDellyPercFreq_90percRecipOverlap,
                                       pmax(SVs$svMantaXPercFreq_90percRecipOverlap,
                                            pmax(SVs$svManta2PercFreq_90percRecipOverlap,
                                                 pmax(SVs$svDellyXPercFreq_90percRecipOverlap,
                                                      SVs$svDelly2PercFreq_90percRecipOverlap, na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T), na.rm=T)/100
  }
  return (SV_MAF)
}

Combine_Tables <- function(table_names) {
  
  res <- data.frame()
  for (i in 1:length(table_names)) {
    res <- rbind(res, read.delim(table_names[[i]], stringsAsFactors = F))
  }
  
  return (res)
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
  }
  
  meta$Status <- ifelse(meta$Relation %in% child_relations, 1, 0)
  return (meta)
}

# Get CH event rate (CH events per bp CNV queried) for individual
Get_SNVs_With_Exonic_Size <- function(CNVs, SNVs, filtered_meta, genes_data, 
                                      with_CNV = T, ds='MSSNG', CNV_freq_filter=1, 
                                      alt_base_filter_threshold=0.9) {
  
  # Determine CNV minor allele frequency
  total_exonic_size <- 0
  return_SNV <- data.frame()
  if (length(CNVs) > 0) {
    CNVs$freq_max <- Determine_CNV_MAF(CNVs)
    
    if (nrow(CNVs) > 0) {
      if (CNVs$Sample.ID[1] %in% filtered_meta$Sample.ID) { # Implicitly apply all filtering
        # Rename startANN and endAnn cols in SVs to match CNVs:
        if (!with_CNV) {
          names(CNVs)[names(CNVs) == 'startAnn'] <- 'STARTAnn'
          names(CNVs)[names(CNVs) == 'endAnn'] <- 'ENDAnn'
        }
        
        CNVs_g <- GRanges(CNVs$chrAnn,
                          IRanges(CNVs$STARTAnn,
                                  CNVs$ENDAnn), "*")
        genes <- genes_data[genes_data$V5 %in% strsplit(paste(CNVs$gene_symbol, collapse = "|"), "\\|")[[1]], ]
        
        if (nrow(genes) > 0) {
          genes_g <- GRanges(genes$V1, IRanges(genes$V2, genes$V3), "*")
          genes_g <- GenomicRanges::reduce(genes_g)
          
          olap <- data.frame(findOverlaps(CNVs_g, genes_g))
          
          if (nrow(olap) > 0) {
            olap$width <- width(pintersect(CNVs_g[olap$queryHits], genes_g[olap$subjectHits]))
            olap <- aggregate(width ~ queryHits, olap, sum)
            CNVs$exonicSize <- 0
            CNVs$exonicSize[olap$queryHits] <- olap$width
            total_exonic_size <- sum(CNVs$exonicSize)
          }
          
          if (nrow(SNVs) > 0) {
            names(SNVs)[11] <- "SampleData"
            SNVs_g <- GRanges(SNVs$CHROM, IRanges(SNVs$start, SNVs$end), "*")
            CNVs_g_redux <- GenomicRanges::reduce(CNVs_g)
            CNVs_g_redux_df <- data.frame(CNVs_g_redux)
            olap <- findOverlaps(SNVs_g, CNVs_g_redux)
            SNVs <- SNVs[unique(olap@from), ]
            SNVs$cnvchrAnn <- CNVs_g_redux_df$seqnames[olap@to]
            SNVs$cnvSTARTAnn <- CNVs_g_redux_df$start[olap@to]
            SNVs$cnvENDAnn <- CNVs_g_redux_df$end[olap@to]
            SNVs$cnvFreq <- CNVs$freq_max[olap@to]
            SNVs$CHfreq <- SNVs$cnvFreq * SNVs$freq_max
            SNVs$exonicSize <- CNVs$exonicSize[olap@to]
            return_SNV <- rbind(return_SNV, SNVs)
          }
        }
      }
    }
  }
  
  return_SNV <- subset(return_SNV, return_SNV$alt_fraction >= alt_base_filter_threshold)
  # Filter ASD associated CNVs and LOFs
  # return_SNV <- Filter_ASD_Associated_CNVs_and_LOFs(return_SNV, ds)
  return (list(return_SNV, total_exonic_size))
}

Process_Parent_Child_CNV_SNV_Data <- function(data_sets, meta_data, genes_data,
                                              with_CNV = T, ds='MSSNG', CNV_freq_filter=1) {
  all_SNVs <- data.frame()
  all_exonic_sizes <- data.frame()
  all_CNVs <- data.frame()
  for (j in 1:length(data_sets)) {
    message(j)
    data <- data_sets[[j]]
    for(i in 1:length(data)){
      childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
      childSNVs <- data.frame(data[[i]]$SNVs)
      
      fatherCNVs <- data.frame(data[[i]]$Father.CNV)
      fatherSNVs <- data.frame(data[[i]]$Father.SNV)
      
      motherCNVs <- data.frame(data[[i]]$Mother.CNV)
      motherSNVs <- data.frame(data[[i]]$Mother.SNV)
      
      if (nrow(childCNVs) > 0) {
        childData <- Get_SNVs_With_Exonic_Size(childCNVs, childSNVs, meta_data, genes_data, with_CNV, ds, CNV_freq_filter)
        childSNVData <- childData[[1]]
        childExonicSize <- data.frame(Sample.ID=childCNVs$Sample.ID[1], exonicSize=childData[[2]])
        all_SNVs <- rbind(all_SNVs, childSNVData)
        all_exonic_sizes <- rbind(all_exonic_sizes, childExonicSize)
        all_CNVs <- rbind(all_CNVs, childCNVs)
      }
      if (nrow(motherCNVs) > 0) {
        motherData <- Get_SNVs_With_Exonic_Size(motherCNVs, motherSNVs, meta_data, genes_data, with_CNV, ds, CNV_freq_filter)
        motherSNVData <- motherData[[1]]
        motherExonicSize <- data.frame(Sample.ID=motherCNVs$Sample.ID[1], exonicSize=motherData[[2]])
        all_SNVs <- rbind(all_SNVs, motherSNVData)
        all_exonic_sizes <- rbind(all_exonic_sizes, motherExonicSize)
        all_CNVs <- rbind(all_CNVs, motherCNVs)
      }
      if (nrow(fatherCNVs) > 0) {
        fatherData <- Get_SNVs_With_Exonic_Size(fatherCNVs, fatherSNVs, meta_data, genes_data, with_CNV, ds, CNV_freq_filter)
        fatherSNVData <- fatherData[[1]]
        fatherExonicSize <- data.frame(Sample.ID=fatherCNVs$Sample.ID[1], exonicSize=fatherData[[2]])
        all_SNVs <- rbind(all_SNVs, fatherSNVData)
        all_exonic_sizes <- rbind(all_exonic_sizes, fatherExonicSize)
        all_CNVs <- rbind(all_CNVs, fatherCNVs)
      }
      
      message(sprintf('%s of %s', i, length(data)))
    }
  }
  
  return (list(all_SNVs, all_exonic_sizes, all_CNVs))
}


gene_exon_data <- data.table::fread("./Data/hg38_refGene_20200708.exon.txt", data.table = F)
### SPARKWGS2 CNV Parent-Proband
SPARKWGS2_metadata_path <- './Data/SPARK_WGS_2_metadata_relfixed.tsv'
SPARKWGS2_meta <- Get_Filtered_Metadata(SPARKWGS2_metadata_path, child_relation = 'proband')

SPARKWGS2_CNV_data <- yaml::yaml.load_file("SPARKWGS2_CH_Data_CNV1P_SNV.yaml")
SPARKWGS2_all_data <- list(SPARKWGS2_CNV_data)
SPARKWGS2_parent_proband_SNVsExonicSizes <- Process_Parent_Child_CNV_SNV_Data(SPARKWGS2_all_data,
                                                                              SPARKWGS2_meta,
                                                                              gene_exon_data,
                                                                              CNV_freq_filter=0.01,
                                                                              ds='SSC')
write_yaml(SPARKWGS2_parent_proband_SNVsExonicSizes, "SPARKWGS2_parent_proband_SNVsExonicSizes.yaml")

