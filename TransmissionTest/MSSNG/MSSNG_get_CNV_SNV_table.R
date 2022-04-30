library(data.table)
library(dplyr)

# setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

Get_Failed_QC_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  # Subjects_to_exclude <- fread('excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  # Subjects_to_exclude <- fread('samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/MSSNG+SSC.CNVs.tsv")
  # ls1 <- fread("MSSNG+SSC.ASD135_LoF.tsv")
  # ls2 <- fread("MSSNG+SSC.CNVs.tsv")
  ls_exclude <- unique(rbind(ls1, ls2))
  
  full_metadata <- fread(file, data.table = F)
  familyID_exclude <- full_metadata$`Family ID`[which(full_metadata$`Sample ID` %in% ls_exclude$Sample)]
  sampleID_exclude <- full_metadata$`Sample ID`[which(full_metadata$`Family ID` %in% familyID_exclude)]
  
  return (sampleID_exclude)
}

Filter_Metadata <- function(file, snvs_IDs) {
  metadata <- read.delim(file, stringsAsFactors = F)
  # Get samples that failed QC
  failed_QC <- Get_Failed_QC_Samples()
  # Get samples with cr CNV
  samples_with_cr_CNV <- Get_cr_CNV_Samples()
  # Get families to exclude that have parents AS the proband as well
  families_to_exclude <- metadata$Family.ID[(metadata$Sample.ID[metadata$Relation == 'proband'] %in% metadata$Mother.ID) |
                                          (metadata$Sample.ID[metadata$Relation == 'proband'] %in% metadata$Father.ID)]
  # Get samples with ASD associated CNVs and LOFs
  ASD_associated_samples <- Get_SampleID_ASD_Associated_CNVs_and_LOFs(file)
  
  # Filter out from metadata
  metadata <- metadata[!metadata$Sample.ID %in% failed_QC, ]
  metadata <- metadata[!metadata$Sample.ID %in% samples_with_cr_CNV, ]
  metadata <- metadata[!metadata$Family.ID %in% families_to_exclude, ]
  metadata <- metadata[!metadata$Sample.ID %in% ASD_associated_samples, ]
  
  # List of samples belonging to families with a proband and both parents data:
  probands <- metadata$Family.ID[which(metadata$Relation %in% c("proband", "affected sibling"))]
  mothers <- metadata$Family.ID[which(metadata$Relation == "mother")]
  fathers <- metadata$Family.ID[which(metadata$Relation == "father")]  
  families <- intersect(intersect(probands, mothers), fathers)  # a vector of Family IDs with a full family dataset
  
  metadata <- metadata[which(metadata$Family.ID %in% families &
                               metadata$Sample.ID %in% snvs_IDs), ]
  
  return (metadata)
}

## Get list of all SNVs:
# path <- "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/TransmissionTest/SNV/MSSNG"
# files <- list.files(path, full.names = T)
# snvs <- data.frame()
# for (i in 1:length(files)) {
#   message(sprintf('Reading snv file: %s of %s', i, length(files)))
#   snvs <- rbind(snvs, fread(files[i], data.table=F))
# }
# write.table(snvs,"MSSNG_SNVs.tsv", sep="\t", row.names=F, quote=F, col.names=T)
snvs <- fread("MSSNG_SNVs.tsv", data.table=F)
message('Done loading snvs')
# snvs <- fread('MSSNG.SNVs.freq.10perc_1.tsv', data.table=F)

## Get filtered metadata file:
metadata <- Filter_Metadata("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/MSSNG_metadata.tsv", unique(snvs$`#Sample`))
# metadata <- Filter_Metadata("MSSNG_metadata.tsv", unique(snvs$`#Sample`))

cnv <- fread("../MSSNG.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("../MSSNG.SVs.freq.1perc.tsv", data.table = F)
# cnv <- fread("MSSNG.CNVs.freq.1perc.tsv", data.table = F)
# sv <- fread("MSSNG.SVs.freq.1perc.tsv", data.table=F)
## Get list of CNVs and SVs that impact CDS:
cnv <- cnv[which(cnv$cds_symbol != "" & cnv$sample %in% metadata$Sample.ID),]
sv <- sv[which(sv$cds_symbol != "" & sv$length <= 10000 
               & sv$sample %in% metadata$Sample.ID),]

#list of samples belonging to families with a proband and both parents with SNV data

## get the list of genes and filter SNVs/Indels data based on the filtered list
cnv.genes <- unique(strsplit(paste(cnv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
sv.genes <- unique(strsplit(paste(sv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
snvs <- snvs[which(snvs$gene_symbol %in% c(cnv.genes, sv.genes)), ]
## Make tables:

Get_CNV_SNV_Table <- function(cnv_and_sv, full_snvs, metadata) {
  table <- data.frame()
  for (i in 1:length(cnv_and_sv)) {
    cnv <- cnv_and_sv[[i]]
    for (cnv_row in 1:nrow(cnv)){
      message(sprintf("%s out of %s", cnv_row, nrow(cnv)))
      cnv.tmp <- cnv[cnv_row, names(cnv)[!names(cnv) %in% names(full_snvs)]]
      inh <- cnv.tmp$Inheritance
      Lof_cds_symbol <- strsplit(cnv.tmp$cds_symbol, "\\|")[[1]] #list of cds_symbols split at |
      
      probandID <- cnv.tmp$sample
      
      ## get non-transmitting and transmitting parent ID
      if (inh == "Paternal") {
        nonInhParentID <- metadata$Mother.ID[which(metadata$Sample.ID==probandID)]
        InhParentID <- metadata$Father.ID[which(metadata$Sample.ID==probandID)]
      } else {  #i.e. if inh == "Maternal"
        nonInhParentID <- metadata$Father.ID[which(metadata$Sample.ID==probandID)]
        InhParentID <- metadata$Mother.ID[which(metadata$Sample.ID==probandID)]
      }
      
      # SNVs that match gene
      snv <- full_snvs[full_snvs$gene_symbol %in% Lof_cds_symbol,]
      
      ## get SNVs of probands, transmitting and non-transmitting parent
      snv.inh.parent <- snv[which(snv$`#Sample` == InhParentID), ]
      snv.noninh.parent <- snv[which(snv$`#Sample` == nonInhParentID), ]
      snv.proband <- snv[which(snv$`#Sample` == probandID), ]
      
      ## get SNPs that found in both parents, as well as homozygous variants in non-inh parent. Then, remove them from the analysis
      snps.to.remove <- c(intersect(snv.inh.parent$"#id", snv.noninh.parent$"#id"),
                          snv.noninh.parent$"#id"[which(snv.noninh.parent$OZYG %in% c("alt-alt", "hom-alt"))])
      snv.noninh.parent <- snv.noninh.parent[which(!snv.noninh.parent$"#id" %in% snps.to.remove), ]
      
      ### retain only SNVs/Indels found in inh parent
      snv.proband <- snv.proband[which(snv.proband$"#id" %in% snv.noninh.parent$"#id"), ]
      snv.tmp <- rbind(snv.noninh.parent, snv.proband)
      
      ## check if snv.tmp is empty
      if (nrow(snv.tmp) == 0) {
        next #move to next cnv row
      }
      snv.tmp[, names(cnv.tmp)] <- cnv.tmp
      table <- bind_rows(table, snv.tmp)
    }
  }
  return (table)
}

table <- Get_CNV_SNV_Table(list(cnv, sv), snvs, metadata)

meta_parentsID <- metadata$Sample.ID[metadata$Relation %in% c('father','mother')]
meta_probandID <- metadata$Sample.ID[metadata$Relation %in% c('proband','affected sibling')]

## add a column to indicate the relation
table$Relation <- NA
table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
table$Relation[which(table$`#Sample` %in% meta_probandID)] <- "proband or affected sibling"
table <- table[which(!table$Relation == "NA"),]
# table <- table[which(!table$CHROM == "chrX"),] #exclude variants on chrX

write.table(table,"MSSNG_CNV_SNV_table.tsv", sep="\t", row.names=F, quote=F, col.names=T)

# plot <- ggplot(data = table, aes(x = effect_priority, fill = Relation))
# plot + geom_bar(stat = "count", position = position_dodge()) +
#   geom_label(stat = "count", aes(label = ..count.., fill = Relation), position = position_dodge(width = 1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
