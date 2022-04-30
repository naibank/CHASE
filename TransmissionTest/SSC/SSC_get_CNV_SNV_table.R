library(data.table)
library(dplyr)

# setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

Get_Failed_QC_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/excluded.qcfailed.samples.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_cr_CNV_Samples <- function() {
  #MSSNG + SSC combined
  Subjects_to_exclude <- fread('/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/samples.with.crCNV.txt', header=F, col.names='Subjects_to_exclude')
  
  return (Subjects_to_exclude$Subjects_to_exclude) # TODO: Fix this, bad design
}

Get_SampleID_ASD_Associated_CNVs_and_LOFs <- function(file) {
  ls1 <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/MSSNG+SSC.ASD135_LoF.tsv")
  ls2 <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/MSSNG+SSC.CNVs.tsv")
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
  US <- metadata$Family.ID[which(metadata$Relation %in% c("unaffected sibling", "other sibling"))]
  mothers <- metadata$Family.ID[which(metadata$Relation == "mother")]
  fathers <- metadata$Family.ID[which(metadata$Relation == "father")]  
  proband_families <- intersect(intersect(probands, mothers), fathers)  # a vector of Family IDs with a full family dataset
  US_families <- intersect(intersect(US, mothers), fathers)  # a vector of Family IDs with a full family dataset
  
  metadata <- metadata[which(metadata$Family.ID %in% c(proband_families,US_families) &
                               metadata$Sample.ID %in% snvs_IDs), ]
  
  return (metadata)
}

## Get list of all SNVs:
# path <- "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/TransmissionTest/SNV/SSC"
# files <- list.files(path, full.names = T)
# snvs <- data.frame()
# for (i in 1:length(files)) {
#   message(sprintf('Reading snv file: %s of %s', i, length(files)))
#   snvs <- rbind(snvs, fread(files[i], data.table=F))
# }
# write.table(snvs,"SSC_SNVs.tsv", sep="\t", row.names=F, quote=F, col.names=T)
snvs <- fread("SSC_SNVs.tsv", data.table=F)
message('Done loading snvs')

## Get filtered metadata file:
metadata <- Filter_Metadata("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/DirectComparison/Data/SSC_metadata.tsv", unique(snvs$`#Sample`))

cnv <- fread("../SSC.CNVs.freq.1perc.tsv", data.table = F)
sv <- fread("../SSC.SVs.freq.1perc.tsv", data.table = F)
## Get list of CNVs and SVs that impact CDS:
cnv <- cnv[which(cnv$cds_symbol != "" & cnv$sample %in% metadata$Sample.ID),]
sv <- sv[which(sv$cds_symbol != "" & sv$length <= 10000 
               & sv$sample %in% metadata$Sample.ID),]

Get_CNV_SNV_Tables <- function(cnv_and_sv, full_snvs, metadata) {
  res <- list()
  for(group in c("proband", "unaffected sibling")){
    ## Make table:
    table <- data.frame()
    
    for (i in 1:length(cnv_and_sv)) {
      cnv.main <- cnv_and_sv[[i]]
      cnv <- cnv.main[cnv.main$sample %in% metadata$Sample.ID[metadata$Relation == group],]
      
      ## get the list of genes and filter SNVs/Indels data based on the filtered list
      cnv.genes <- unique(strsplit(paste(cnv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
      snvs <- full_snvs[which(full_snvs$gene_symbol %in% cnv.genes), ]
      
      for (cnv_row in 1:nrow(cnv)){
        message(sprintf("%s out of %s", cnv_row, nrow(cnv)))
        cnv.tmp <- cnv[cnv_row, names(cnv)[!names(cnv) %in% names(snvs)]]
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
        
        snv <- snvs[which(snvs$gene_symbol %in% Lof_cds_symbol),]
        
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
    
    message(sprintf('done %s', group))
    res[[group]] <- table
  }
  return (res)
}

## remove CRVs from tables
tables <- Get_CNV_SNV_Tables(list(cnv,sv), snvs, metadata)
proband_table <- tables[[1]]
unaffectedsib_table <- tables[[2]]

meta_parentsID <- metadata$Sample.ID[metadata$Relation %in% c('father','mother')]
meta_probandID <- metadata$Sample.ID[metadata$Relation %in% c('proband','affected sibling')]
meta_USID <- metadata$Sample.ID[metadata$Relation %in% c('unaffected sibling','other sibling')]

## add a column to indicate the relation
proband_table$Relation <- NA
proband_table$Relation[which(proband_table$`#Sample` %in% meta_parentsID)] <- "parent"
proband_table$Relation[which(proband_table$`#Sample` %in% meta_probandID)] <- "proband or affected sibling"
proband_table <- proband_table[which(!proband_table$Relation == "NA"),]

unaffectedsib_table$Relation <- NA
unaffectedsib_table$Relation[which(unaffectedsib_table$`#Sample` %in% meta_parentsID)] <- "parent"
unaffectedsib_table$Relation[which(unaffectedsib_table$`#Sample` %in% meta_USID)] <- "unaffected sibling"

write.table(proband_table, "SSC_proband_CNV_SNV_table.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
write.table(unaffectedsib_table, "SSC_unaffected_sibling_CNV_SNV_table.tsv",  sep="\t", row.names=F, quote=F, col.names=T)
