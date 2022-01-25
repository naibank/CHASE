#Required Packages
#BiocManager::install("GenomicRanges")
#BiocManager::install("Repitools")

library(GenomicRanges)
library(Repitools)
library(yaml)
library(data.table)

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

Get_SNVs_Parent_Comparison <- function(file_path, CNV_freq_filter=1, SNV_freq_filter=1,
                                       alt_base_filter_threshold = 0.9) {
  data <- yaml::yaml.load_file(file_path)
  
  childSNVs.in.maternalCNVs <- data.frame()
  childSNVs.in.paternalCNVs <- data.frame()
  motherSNVs.in.maternalCNVs <- data.frame()
  fatherSNVs.in.paternalCNVs <- data.frame()

  maternalCount = 0
  paternalCount = 0
  
  motherCount = 0
  fatherCount = 0
  childBothPaternalMaternalCount = 0
  
  for(i in 1:length(data)){
    childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
    childSNVs <- data.frame(data[[i]]$SNVs)
    
    fatherCNVs <- data.frame(data[[i]]$Father.CNV)
    fatherSNVs <- data.frame(data[[i]]$Father.SNV)
    
    motherCNVs <- data.frame(data[[i]]$Mother.CNV)
    motherSNVs <- data.frame(data[[i]]$Mother.SNV)
    
    #Determine CNV minor allele frequency
    if (length(childCNVs) > 0) {
      childCNVs$cnvFreq <- Determine_CNV_MAF(childCNVs)
      childCNVs <-  childCNVs[which(childCNVs$cnvFreq < CNV_freq_filter), ]
    }
    if (length(fatherCNVs) > 0) {
      fatherCNVs$cnvFreq <- Determine_CNV_MAF(fatherCNVs)
      fatherCNVs <-  fatherCNVs[which(fatherCNVs$cnvFreq < CNV_freq_filter), ]
    }
    if (length(motherCNVs) > 0) {
      motherCNVs$cnvFreq <- Determine_CNV_MAF(motherCNVs)
      motherCNVs <-  motherCNVs[which(motherCNVs$cnvFreq < CNV_freq_filter), ]
    }
    
    childSNVs <- childSNVs[which(childSNVs$freq_max < SNV_freq_filter), ]
    fatherSNVs <- fatherSNVs[which(fatherSNVs$freq_max < SNV_freq_filter), ]
    motherSNVs <- motherSNVs[which(motherSNVs$freq_max < SNV_freq_filter), ]
    
    motherCount = motherCount + as.numeric("Maternal" %in% childCNVs$Inheritance)
    fatherCount = fatherCount + as.numeric("Paternal" %in% childCNVs$Inheritance)
    childBothPaternalMaternalCount = childBothPaternalMaternalCount + as.numeric("Paternal" %in% childCNVs$Inheritance &
                                                                                       "Maternal" %in% childCNVs$Inheritance)
    if(!("Paternal" %in% childCNVs$Inheritance & "Maternal" %in% childCNVs$Inheritance)){
      if(sum(childCNVs$Inheritance == "Maternal") > 0){
        maternalCount = maternalCount + 1 
        
        message(sprintf("%s %s %s", i, motherCount-childBothPaternalMaternalCount, maternalCount))
        child.motherCNVs <- GRanges(childCNVs$chrAnn[childCNVs$Inheritance == "Maternal"],
                                      IRanges(childCNVs$STARTAnn[childCNVs$Inheritance == "Maternal"],
                                              childCNVs$ENDAnn[childCNVs$Inheritance == "Maternal"]), "*")
        
        if(nrow(childSNVs) > 0){
          names(childSNVs)[11] <- "SampleData"
          
          childSNVs.g <- GRanges(childSNVs$CHROM, IRanges(childSNVs$start, childSNVs$end), "*")
          olap <- findOverlaps(childSNVs.g, child.motherCNVs)
          childSNVs <- childSNVs[unique(olap@from), ]
          childSNVs$cnvFreq <- childCNVs$cnvFreq[olap@to]
          childSNVs$CHfreq <- childSNVs$cnvFreq * childSNVs$freq_max
          childSNVs.in.maternalCNVs <- rbind(childSNVs.in.maternalCNVs, childSNVs)
        }
        
        if(nrow(motherCNVs) > 0 & nrow(motherSNVs) > 0){
          names(motherSNVs)[11] <- "SampleData"
          motherCNVs.g <- GRanges(motherCNVs$chrAnn,
                                IRanges(motherCNVs$STARTAnn,
                                        motherCNVs$ENDAnn), "*")
          motherSNVs <- motherSNVs[,names(motherSNVs)[names(motherSNVs) != "Inheritance.SNV"]]
          motherSNVs.g <- GRanges(motherSNVs$CHROM, IRanges(motherSNVs$start, motherSNVs$end), "*")
          olap <- findOverlaps(motherSNVs.g, motherCNVs.g)
          motherSNVs <- motherSNVs[unique(olap@from), ]
          motherSNVs$cnvFreq <- motherCNVs$cnvFreq[olap@to]
          motherSNVs$CHfreq <- motherSNVs$cnvFreq * motherSNVs$freq_max
          motherSNVs.in.maternalCNVs <- rbind(motherSNVs.in.maternalCNVs, motherSNVs)
        }
      }
      
      if(sum(childCNVs$Inheritance == "Paternal") > 0){ 
        paternalCount = paternalCount + 1
        child.fatherCNVs <- GRanges(childCNVs$chrAnn[childCNVs$Inheritance == "Paternal"],
                                      IRanges(childCNVs$STARTAnn[childCNVs$Inheritance == "Paternal"],
                                              childCNVs$ENDAnn[childCNVs$Inheritance == "Paternal"]), "*")

        if(nrow(childSNVs) > 0){
          names(childSNVs)[11] <- "SampleData"
          
          childSNVs.g <- GRanges(childSNVs$CHROM, IRanges(childSNVs$start, childSNVs$end), "*")
          olap <- findOverlaps(childSNVs.g, child.fatherCNVs)
          childSNVs <- childSNVs[unique(olap@from), ]
          childSNVs$cnvFreq <- childCNVs$cnvFreq[olap@to]
          childSNVs$CHfreq <- childSNVs$cnvFreq * childSNVs$freq_max
          childSNVs.in.paternalCNVs <- rbind(childSNVs.in.paternalCNVs, childSNVs)
        }
        
        if(nrow(fatherCNVs) > 0 & nrow(fatherSNVs) > 0){
          names(fatherSNVs)[11] <- "SampleData"
          fatherCNVs.g <- GRanges(fatherCNVs$chrAnn,
                                IRanges(fatherCNVs$STARTAnn,
                                        fatherCNVs$ENDAnn), "*")
          fatherSNVs <- fatherSNVs[,names(fatherSNVs)[names(fatherSNVs) != "Inheritance.SNV"]]
          
          fatherSNVs.g <- GRanges(fatherSNVs$CHROM, IRanges(fatherSNVs$start, fatherSNVs$end), "*")
          olap <- findOverlaps(fatherSNVs.g, fatherCNVs.g)
          fatherSNVs <- fatherSNVs[unique(olap@from), ]
          fatherSNVs$cnvFreq <- fatherCNVs$cnvFreq[olap@to]
          fatherSNVs$CHfreq <- fatherSNVs$cnvFreq * fatherSNVs$freq_max
          fatherSNVs.in.paternalCNVs <- rbind(fatherSNVs.in.paternalCNVs, fatherSNVs)
        }
      }
    }
    message(i)
  }
  
  childSNVs.in.paternalCNVs <- subset(childSNVs.in.paternalCNVs, childSNVs.in.paternalCNVs$alt_fraction >= alt_base_filter_threshold)
  childSNVs.in.maternalCNVs <- subset(childSNVs.in.maternalCNVs, childSNVs.in.maternalCNVs$alt_fraction >= alt_base_filter_threshold)
  fatherSNVs.in.paternalCNVs <- subset(fatherSNVs.in.paternalCNVs, fatherSNVs.in.paternalCNVs$alt_fraction >= alt_base_filter_threshold)
  motherSNVs.in.maternalCNVs <- subset(motherSNVs.in.maternalCNVs, motherSNVs.in.maternalCNVs$alt_fraction >= alt_base_filter_threshold)
  
  return (list(childSNVs.in.paternalCNVs, 
               childSNVs.in.maternalCNVs, 
               fatherSNVs.in.paternalCNVs, 
               motherSNVs.in.maternalCNVs,
               paternalCount,
               maternalCount))
}

Get_SNVs_Sibling_Comparison <- function(file_path, CNV_freq_filter=1, SNV_freq_filter=1, 
                                        alt_base_filter_threshold = 0.9) {
  data <- yaml::yaml.load_file(file_path)
  
  childSNVs_combined <- data.frame()
  SNVs_to_add <- data.frame()
  
  total_exonic_size = 0
  genes_data <- data.table::fread("hg38_refGene_20200708.exon.txt", data.table = F)
  
  for(i in 1:length(data)){
    childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
    childSNVs <- data.frame(data[[i]]$SNVs)
    
    if(length(childCNVs) > 0) {
      childCNVs$cnvFreq <- Determine_CNV_MAF(childCNVs)
      childCNVs <-  childCNVs[which(childCNVs$cnvFreq < CNV_freq_filter), ]
    }

    childSNVs <- childSNVs[which(childSNVs$freq_max < SNV_freq_filter), ]
    
    if (nrow(childCNVs) > 0) {
      
      childCNVs_g <- GRanges(childCNVs$chrAnn,
                             IRanges(childCNVs$STARTAnn,
                                     childCNVs$ENDAnn), "*")
      
      genes <- genes_data[genes_data$V5 %in% strsplit(paste(childCNVs$gene_symbol, collapse = "|"), "\\|")[[1]], ]
      
      genes.g <- GRanges(genes$V1, IRanges(genes$V2, genes$V3), "*")
      genes.g <- GenomicRanges::reduce(genes.g)
      
      olap <- data.frame(findOverlaps(childCNVs_g, genes.g))
      olap$width <- width(pintersect(childCNVs_g[olap$queryHits], genes.g[olap$subjectHits]))
      olap <- aggregate(width ~ queryHits, olap, sum)
      childCNVs$exonicSize <- 0
      childCNVs$exonicSize[olap$queryHits] <- olap$width
      total_exonic_size <- sum(childCNVs$exonicSize, total_exonic_size)
      
      if(nrow(childSNVs) > 0){
        names(childSNVs)[11] <- "SampleData"
        
        childSNVs_g <- GRanges(childSNVs$CHROM, IRanges(childSNVs$start, childSNVs$end), "*")
        olap <- findOverlaps(childSNVs_g, childCNVs_g)
        childSNVs <- childSNVs[unique(olap@from), ]
        childSNVs$cnvFreq <- childCNVs$cnvFreq[olap@to]
        childSNVs$CHfreq <- childSNVs$cnvFreq * childSNVs$freq_max
        childSNVs$cnvInheritance <- childCNVs$Inheritance[olap@to]
        childSNVs_combined <- rbind(childSNVs_combined, childSNVs)
      }
    }
    message(i)
  }
  
  # Exclude SNVs inherited from the deletion transmitting parent via alt_fraction column
  childSNVs_combined <- subset(childSNVs_combined, 
                               childSNVs_combined$alt_fraction >= alt_base_filter_threshold)
  
  
  return (list(childSNVs_combined, total_exonic_size))
}

##### ILMN DATA ##### (1111 CH events)

ILMN_SNVs_data <- Get_SNVs_Parent_Comparison("MSSNG_ILMN_CH_Data_CNV10P_SNV.yaml", CNV_freq_filter=0.01, SNV_freq_filter=1)

ILMN_probandSNVs_in_paternalCNVs <- ILMN_SNVs_data[[1]]
ILMN_probandSNVs_in_maternalCNVs <- ILMN_SNVs_data[[2]]
ILMN_fatherSNVs_in_paternalCNVs <- ILMN_SNVs_data[[3]]
ILMN_motherSNVs_in_maternalCNVs <- ILMN_SNVs_data[[4]]
ILMN_paternalCNV_inh_count <- ILMN_SNVs_data[[5]]
ILMN_maternalCNV_inh_count <- ILMN_SNVs_data[[6]]

# Write SNVs to tsv
write.table(ILMN_probandSNVs_in_paternalCNVs, "../DT/ILMN_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(ILMN_probandSNVs_in_maternalCNVs, "../DT/ILMN_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(ILMN_fatherSNVs_in_paternalCNVs, "../DT/ILMN_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(ILMN_motherSNVs_in_maternalCNVs, "../DT/ILMN_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

##### CG DATA ##### (120 CH events)

CG_SNVs_data <- Get_SNVs_Parent_Comparison("MSSNG_CG_CH_Data_CNV10P_SNV.yaml", CNV_freq_filter=0.01, SNV_freq_filter=1)

CG_probandSNVs_in_paternalCNVs <- CG_SNVs_data[[1]]
CG_probandSNVs_in_paternalCNVs[which(CG_probandSNVs_in_paternalCNVs$X.Sample == "-717640"), "X.Sample"] <- "5-0003-003" # Sample w error in ID 
CG_probandSNVs_in_maternalCNVs <- CG_SNVs_data[[2]]
CG_fatherSNVs_in_paternalCNVs <- CG_SNVs_data[[3]]
CG_motherSNVs_in_maternalCNVs <- CG_SNVs_data[[4]]
CG_paternalCNV_inh_count <- CG_SNVs_data[[5]]
CG_maternalCNV_inh_count <- CG_SNVs_data[[6]]

# Write SNVs to tsv
write.table(CG_probandSNVs_in_paternalCNVs, "../DT/CG_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(CG_probandSNVs_in_maternalCNVs, "../DT/CG_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(CG_fatherSNVs_in_paternalCNVs, "../DT/CG_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(CG_motherSNVs_in_maternalCNVs, "../DT/CG_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

##### SSC DATA ##### (1510 CH events)

SSC_SNVs_data <- Get_SNVs_Parent_Comparison("SSC_CH_Data_CNV10P_SNV.yaml", CNV_freq_filter=0.01, SNV_freq_filter=1)

SSC_probandSNVs_in_paternalCNVs <- SSC_SNVs_data[[1]]
SSC_probandSNVs_in_maternalCNVs <- SSC_SNVs_data[[2]]
SSC_fatherSNVs_in_paternalCNVs <- SSC_SNVs_data[[3]]
SSC_motherSNVs_in_maternalCNVs <- SSC_SNVs_data[[4]]
SSC_paternalCNV_inh_count <- SSC_SNVs_data[[5]]
SSC_maternalCNV_inh_count <- SSC_SNVs_data[[6]]

# Write SNVs to tsv
write.table(SSC_probandSNVs_in_paternalCNVs, "../DT/SSC_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_probandSNVs_in_maternalCNVs, "../DT/SSC_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_fatherSNVs_in_paternalCNVs, "../DT/SSC_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_motherSNVs_in_maternalCNVs, "../DT/SSC_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

##### SSC DATA - Proband vs. Unaffected Siblings ##### (X CH events)

SSC_probandUN_comparison_results <- Get_SNVs_Sibling_Comparison("SSC_CH_Data_CNV10P_SNV.yaml", CNV_freq_filter=0.01, SNV_freq_filter=1)
SSC_probandUNSNVs <- SSC_probandUN_comparison_results[[1]]
SSC_probandUN_exonic_size <- SSC_probandUN_comparison_results[[2]]

SSC_unaffectedSiblings_comparison_results <- Get_SNVs_Sibling_Comparison("SSC_CH.unaffectedSiblings_Data_CNV10P_SNV.yaml", CNV_freq_filter=0.01, SNV_freq_filter=1)
SSC_unaffectedSiblingsSNVs <- SSC_unaffectedSiblings_comparison_results[[1]]
SSC_unaffectedSiblings_exonic_size <- SSC_unaffectedSiblings_comparison_results[[2]]

# Write SNVs to tsv
write.table(SSC_probandUNSNVs, "../DT/SSC_ProbandUNSNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_unaffectedSiblingsSNVs, "../DT/SSC_UnaffectedSiblingsSNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

# Optional: Restrict by CH event frequency of found SNVs; off when running script all at once
restrict <- F
if (restrict) {
  CH_Freq_Threshold <- 0.001
  
  ILMN_probandSNVs_in_paternalCNVs <- subset(ILMN_probandSNVs_in_paternalCNVs, ILMN_probandSNVs_in_paternalCNVs$CHfreq <= CH_Freq_Threshold)
  ILMN_probandSNVs_in_maternalCNVs <- subset(ILMN_probandSNVs_in_maternalCNVs, ILMN_probandSNVs_in_maternalCNVs$CHfreq <= CH_Freq_Threshold)
  ILMN_fatherSNVs_in_paternalCNVs <- subset(ILMN_fatherSNVs_in_paternalCNVs, ILMN_fatherSNVs_in_paternalCNVs$CHfreq <= CH_Freq_Threshold)
  ILMN_motherSNVs_in_maternalCNVs <- subset(ILMN_motherSNVs_in_maternalCNVs, ILMN_motherSNVs_in_maternalCNVs$CHfreq <= CH_Freq_Threshold)
  
  CG_probandSNVs_in_paternalCNVs <- subset(CG_probandSNVs_in_paternalCNVs, CG_probandSNVs_in_paternalCNVs$CHfreq <= CH_Freq_Threshold)
  CG_probandSNVs_in_maternalCNVs <- subset(CG_probandSNVs_in_maternalCNVs, CG_probandSNVs_in_maternalCNVs$CHfreq <= CH_Freq_Threshold)
  CG_fatherSNVs_in_paternalCNVs <- subset(CG_fatherSNVs_in_paternalCNVs, CG_fatherSNVs_in_paternalCNVs$CHfreq <= CH_Freq_Threshold)
  CG_motherSNVs_in_maternalCNVs <- subset(CG_motherSNVs_in_maternalCNVs, CG_motherSNVs_in_maternalCNVs$CHfreq <= CH_Freq_Threshold)
  
  SSC_probandSNVs_in_paternalCNVs <- subset(SSC_probandSNVs_in_paternalCNVs, SSC_probandSNVs_in_paternalCNVs$CHfreq <= CH_Freq_Threshold)
  SSC_probandSNVs_in_maternalCNVs <- subset(SSC_probandSNVs_in_maternalCNVs, SSC_probandSNVs_in_maternalCNVs$CHfreq <= CH_Freq_Threshold)
  SSC_fatherSNVs_in_paternalCNVs <- subset(SSC_fatherSNVs_in_paternalCNVs, SSC_fatherSNVs_in_paternalCNVs$CHfreq <= CH_Freq_Threshold)
  SSC_motherSNVs_in_maternalCNVs <- subset(SSC_motherSNVs_in_maternalCNVs, SSC_motherSNVs_in_maternalCNVs$CHfreq <= CH_Freq_Threshold)
  
  SSC_probandUNSNVs <- subset(SSC_probandUNSNVs, SSC_probandUNSNVs$CHfreq <= CH_Freq_Threshold)
  SSC_unaffectedSiblingsSNVs <- subset(SSC_unaffectedSiblingsSNVs, SSC_unaffectedSiblingsSNVs$CHfreq <= CH_Freq_Threshold)
}
###########################################################################################################################################################

# SNV COUNT BY VARIANT TYPE

Get_SNV_Count_By_Var_Type <- function(file_path, loaded_data=NULL) {
  # Read SNVs tsv
  if (is.null(loaded_data)) {
    SNVs_data <- read.delim(file_path, stringsAsFactors = F)
  }
  else {
    SNVs_data <- loaded_data
  }

  # All variants
  SNVs_all <- length(SNVs_data$X.Sample[!duplicated(paste(SNVs_data$X.Sample, SNVs_data$X.id))])                
  SNVs_all_count <- length(unique(SNVs_data$X.Sample))
  
  # Synonymous SNVs
  SNVs_syn <- nrow(unique(SNVs_data[which(SNVs_data$effect_priority == "synonymous SNV"), c("X.Sample", "X.id")]))           
  SNVs_syn_count <- length(unique(SNVs_data$X.Sample[which(SNVs_data$effect_priority == "synonymous SNV")]))
  
  # Missense SNVs
  SNVs_mis <- nrow(unique(SNVs_data[which(SNVs_data$effect_priority == "nonsynonymous SNV"), c("X.Sample", "X.id")]))           
  SNVs_mis_count <- length(unique(SNVs_data$X.Sample[which(SNVs_data$effect_priority == "nonsynonymous SNV")]))
  
  # Loss of Function SNVs
  SNVs_lof <- nrow(unique(SNVs_data[which(SNVs_data$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                            SNVs_data$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
  SNVs_lof_count <- length(unique(SNVs_data$X.Sample[which(SNVs_data$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                             SNVs_data$typeseq_priority == "splicing")]))
  
  return (list(list(all=SNVs_all, syn=SNVs_syn, mis=SNVs_mis, lof=SNVs_lof), 
          list(all=SNVs_all_count, syn=SNVs_syn_count, mis=SNVs_mis_count, lof=SNVs_lof_count)))
}

##### SNV Counts - Stratify by proband sex ######
Get_IDs_By_Sex <- function(metadata_file) {
  metadata <- as.data.frame(fread(metadata_file))
  males <- metadata[which(metadata$Sex == 'male' & metadata$Relation == 'proband'), c("Sample ID", "Father ID", "Mother ID")]
  male_IDs <- males$'Sample ID'
  father_of_male_IDs <- males$'Father ID'
  mother_of_male_IDs <- males$'Mother ID'
  
  females <- metadata[which(metadata$Sex == 'female' & metadata$Relation == 'proband'), c("Sample ID", "Father ID", "Mother ID")]
  female_IDs <- females$'Sample ID'
  father_of_female_IDs <- females$'Father ID'
  mother_of_female_IDs <- females$'Mother ID'
  
  return (list(male_IDs, father_of_male_IDs, mother_of_male_IDs, 
               female_IDs, father_of_female_IDs, mother_of_female_IDs))
}

##### MSSNG Data Counts ######
MSSNG_sex_IDs <- Get_IDs_By_Sex("MSSNG_metadata.tsv")
MSSNG_male_IDs <- MSSNG_sex_IDs[[1]]
MSSNG_father_of_male_IDs <- MSSNG_sex_IDs[[2]]
MSSNG_mother_of_male_IDs <- MSSNG_sex_IDs[[3]]

MSSNG_female_IDs <- MSSNG_sex_IDs[[4]]
MSSNG_father_of_female_IDs <- MSSNG_sex_IDs[[5]]
MSSNG_mother_of_female_IDs <- MSSNG_sex_IDs[[6]]

ILMN_male_probandSNVs_in_paternal_CNVs <- ILMN_probandSNVs_in_paternalCNVs[which(ILMN_probandSNVs_in_paternalCNVs$X.Sample %in% MSSNG_male_IDs), ]
ILMN_male_probandSNVs_in_maternal_CNVs <- ILMN_probandSNVs_in_maternalCNVs[which(ILMN_probandSNVs_in_maternalCNVs$X.Sample %in% MSSNG_male_IDs), ]
ILMN_father_of_male_SNVs_in_paternal_CNVs <- ILMN_fatherSNVs_in_paternalCNVs[which(ILMN_fatherSNVs_in_paternalCNVs$X.Sample %in% MSSNG_father_of_male_IDs), ]
ILMN_mother_of_male_SNVs_in_maternal_CNVs <- ILMN_motherSNVs_in_maternalCNVs[which(ILMN_motherSNVs_in_maternalCNVs$X.Sample %in% MSSNG_mother_of_male_IDs), ]

ILMN_female_probandSNVs_in_paternal_CNVs <- ILMN_probandSNVs_in_paternalCNVs[which(ILMN_probandSNVs_in_paternalCNVs$X.Sample %in% MSSNG_female_IDs), ]
ILMN_female_probandSNVs_in_maternal_CNVs <- ILMN_probandSNVs_in_maternalCNVs[which(ILMN_probandSNVs_in_maternalCNVs$X.Sample %in% MSSNG_female_IDs), ]
ILMN_father_of_female_SNVs_in_paternal_CNVs <- ILMN_fatherSNVs_in_paternalCNVs[which(ILMN_fatherSNVs_in_paternalCNVs$X.Sample %in% MSSNG_father_of_female_IDs), ]
ILMN_mother_of_female_SNVs_in_maternal_CNVs <- ILMN_motherSNVs_in_maternalCNVs[which(ILMN_motherSNVs_in_maternalCNVs$X.Sample %in% MSSNG_mother_of_female_IDs), ]

CG_male_probandSNVs_in_paternal_CNVs <- CG_probandSNVs_in_paternalCNVs[which(CG_probandSNVs_in_paternalCNVs$X.Sample %in% MSSNG_male_IDs), ]
CG_male_probandSNVs_in_maternal_CNVs <- CG_probandSNVs_in_maternalCNVs[which(CG_probandSNVs_in_maternalCNVs$X.Sample %in% MSSNG_male_IDs), ]
CG_father_of_male_SNVs_in_paternal_CNVs <- CG_fatherSNVs_in_paternalCNVs[which(CG_fatherSNVs_in_paternalCNVs$X.Sample %in% MSSNG_father_of_male_IDs), ]
CG_mother_of_male_SNVs_in_maternal_CNVs <- CG_motherSNVs_in_maternalCNVs[which(CG_motherSNVs_in_maternalCNVs$X.Sample %in% MSSNG_mother_of_male_IDs), ]

CG_female_probandSNVs_in_paternal_CNVs <- CG_probandSNVs_in_paternalCNVs[which(CG_probandSNVs_in_paternalCNVs$X.Sample %in% MSSNG_female_IDs), ]
CG_female_probandSNVs_in_maternal_CNVs <- CG_probandSNVs_in_maternalCNVs[which(CG_probandSNVs_in_maternalCNVs$X.Sample %in% MSSNG_female_IDs), ]
CG_father_of_female_SNVs_in_paternal_CNVs <- CG_fatherSNVs_in_paternalCNVs[which(CG_fatherSNVs_in_paternalCNVs$X.Sample %in% MSSNG_father_of_female_IDs), ]
CG_mother_of_female_SNVs_in_maternal_CNVs <- CG_motherSNVs_in_maternalCNVs[which(CG_motherSNVs_in_maternalCNVs$X.Sample %in% MSSNG_mother_of_female_IDs), ]

ILMN_male_paternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_male_probandSNVs_in_paternal_CNVs)
ILMN_male_paternal_SNV_counts <- ILMN_male_paternal_both_counts[[1]]
ILMN_male_paternal_counts <- ILMN_male_paternal_both_counts[[2]]

ILMN_male_maternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_male_probandSNVs_in_maternal_CNVs)
ILMN_male_maternal_SNV_counts <- ILMN_male_maternal_both_counts[[1]]
ILMN_male_maternal_counts <- ILMN_male_maternal_both_counts[[2]]

ILMN_father_of_male_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_father_of_male_SNVs_in_paternal_CNVs)
ILMN_father_of_male_SNV_counts <- ILMN_father_of_male_both_counts[[1]]
ILMN_father_of_male_counts <- ILMN_father_of_male_both_counts[[2]]

ILMN_mother_of_male_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_mother_of_male_SNVs_in_maternal_CNVs)
ILMN_mother_of_male_SNV_counts <- ILMN_mother_of_male_both_counts[[1]]
ILMN_mother_of_male_counts <- ILMN_mother_of_male_both_counts[[2]]

ILMN_female_paternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_female_probandSNVs_in_paternal_CNVs)
ILMN_female_paternal_SNV_counts <- ILMN_female_paternal_both_counts[[1]]
ILMN_female_paternal_counts <- ILMN_female_paternal_both_counts[[2]]

ILMN_female_maternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_female_probandSNVs_in_maternal_CNVs)
ILMN_female_maternal_SNV_counts <- ILMN_female_maternal_both_counts[[1]]
ILMN_female_maternal_counts <- ILMN_female_maternal_both_counts[[2]]

ILMN_father_of_female_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_father_of_female_SNVs_in_paternal_CNVs)
ILMN_father_of_female_SNV_counts <- ILMN_father_of_female_both_counts[[1]]
ILMN_father_of_female_counts <- ILMN_father_of_female_both_counts[[2]]

ILMN_mother_of_female_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = ILMN_mother_of_female_SNVs_in_maternal_CNVs)
ILMN_mother_of_female_SNV_counts <- ILMN_mother_of_female_both_counts[[1]]
ILMN_mother_of_female_counts <- ILMN_mother_of_female_both_counts[[2]]

ILMN_fullSNV_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                male_proband_paternal_del = c(ILMN_male_paternal_SNV_counts$all, ILMN_male_paternal_SNV_counts$syn, ILMN_male_paternal_SNV_counts$mis, ILMN_male_paternal_SNV_counts$lof),
                                father_of_male = c(ILMN_father_of_male_SNV_counts$all, ILMN_father_of_male_SNV_counts$syn, ILMN_father_of_male_SNV_counts$mis, ILMN_father_of_male_SNV_counts$lof),
                                female_proband_paternal_del = c(ILMN_female_paternal_SNV_counts$all, ILMN_female_paternal_SNV_counts$syn, ILMN_female_paternal_SNV_counts$mis, ILMN_female_paternal_SNV_counts$lof),
                                father_of_female = c(ILMN_father_of_female_SNV_counts$all, ILMN_father_of_female_SNV_counts$syn, ILMN_father_of_female_SNV_counts$mis, ILMN_father_of_female_SNV_counts$lof),
                                male_proband_maternal_del = c(ILMN_male_maternal_SNV_counts$all, ILMN_male_maternal_SNV_counts$syn, ILMN_male_maternal_SNV_counts$mis, ILMN_male_maternal_SNV_counts$lof),
                                mother_of_male = c(ILMN_mother_of_male_SNV_counts$all, ILMN_mother_of_male_SNV_counts$syn, ILMN_mother_of_male_SNV_counts$mis, ILMN_mother_of_male_SNV_counts$lof),
                                female_proband_maternal_del = c(ILMN_female_maternal_SNV_counts$all, ILMN_female_maternal_SNV_counts$syn, ILMN_female_maternal_SNV_counts$mis, ILMN_female_maternal_SNV_counts$lof),
                                mother_of_female = c(ILMN_mother_of_female_SNV_counts$all, ILMN_mother_of_female_SNV_counts$syn, ILMN_mother_of_female_SNV_counts$mis, ILMN_mother_of_female_SNV_counts$lof))

ILMN_fullSNV_comparison_df <- as.data.frame(ILMN_fullSNV_comparison)

ILMN_full_indiv_count_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                         male_proband_paternal_del = c(ILMN_male_paternal_counts$all, ILMN_male_paternal_counts$syn, ILMN_male_paternal_counts$mis, ILMN_male_paternal_counts$lof),
                                         father_of_male = c(ILMN_father_of_male_counts$all, ILMN_father_of_male_counts$syn, ILMN_father_of_male_counts$mis, ILMN_father_of_male_counts$lof),
                                         female_proband_paternal_del = c(ILMN_female_paternal_counts$all, ILMN_female_paternal_counts$syn, ILMN_female_paternal_counts$mis, ILMN_female_paternal_counts$lof),
                                         father_of_female = c(ILMN_father_of_female_counts$all, ILMN_father_of_female_counts$syn, ILMN_father_of_female_counts$mis, ILMN_father_of_female_counts$lof),
                                         male_proband_maternal_del = c(ILMN_male_maternal_counts$all, ILMN_male_maternal_counts$syn, ILMN_male_maternal_counts$mis, ILMN_male_maternal_counts$lof),
                                         mother_of_male = c(ILMN_mother_of_male_counts$all, ILMN_mother_of_male_counts$syn, ILMN_mother_of_male_counts$mis, ILMN_mother_of_male_counts$lof),
                                         female_proband_maternal_del = c(ILMN_female_maternal_counts$all, ILMN_female_maternal_counts$syn, ILMN_female_maternal_counts$mis, ILMN_female_maternal_counts$lof),
                                         mother_of_female = c(ILMN_mother_of_female_counts$all, ILMN_mother_of_female_counts$syn, ILMN_mother_of_female_counts$mis, ILMN_mother_of_female_counts$lof))

ILMN_full_indiv_count_comparison_df <- as.data.frame(ILMN_full_indiv_count_comparison)

CG_male_paternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_male_probandSNVs_in_paternal_CNVs)
CG_male_paternal_SNV_counts <- CG_male_paternal_both_counts[[1]]
CG_male_paternal_counts <- CG_male_paternal_both_counts[[2]]

CG_male_maternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_male_probandSNVs_in_maternal_CNVs)
CG_male_maternal_SNV_counts <- CG_male_maternal_both_counts[[1]]
CG_male_maternal_counts <- CG_male_maternal_both_counts[[2]]

CG_father_of_male_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_father_of_male_SNVs_in_paternal_CNVs)
CG_father_of_male_SNV_counts <- CG_father_of_male_both_counts[[1]]
CG_father_of_male_counts <- CG_father_of_male_both_counts[[2]]

CG_mother_of_male_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_mother_of_male_SNVs_in_maternal_CNVs)
CG_mother_of_male_SNV_counts <- CG_mother_of_male_both_counts[[1]]
CG_mother_of_male_counts <- CG_mother_of_male_both_counts[[2]]

CG_female_paternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_female_probandSNVs_in_paternal_CNVs)
CG_female_paternal_SNV_counts <- CG_female_paternal_both_counts[[1]]
CG_female_paternal_counts <- CG_female_paternal_both_counts[[2]]

CG_female_maternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_female_probandSNVs_in_maternal_CNVs)
CG_female_maternal_SNV_counts <- CG_female_maternal_both_counts[[1]]
CG_female_maternal_counts <- CG_female_maternal_both_counts[[2]]

CG_father_of_female_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_father_of_female_SNVs_in_paternal_CNVs)
CG_father_of_female_SNV_counts <- CG_father_of_female_both_counts[[1]]
CG_father_of_female_counts <- CG_father_of_female_both_counts[[2]]

CG_mother_of_female_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = CG_mother_of_female_SNVs_in_maternal_CNVs)
CG_mother_of_female_SNV_counts <- CG_mother_of_female_both_counts[[1]]
CG_mother_of_female_counts <- CG_mother_of_female_both_counts[[2]]

CG_fullSNV_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                              male_proband_paternal_del = c(CG_male_paternal_SNV_counts$all, CG_male_paternal_SNV_counts$syn, CG_male_paternal_SNV_counts$mis, CG_male_paternal_SNV_counts$lof),
                              father_of_male = c(CG_father_of_male_SNV_counts$all, CG_father_of_male_SNV_counts$syn, CG_father_of_male_SNV_counts$mis, CG_father_of_male_SNV_counts$lof),
                              female_proband_paternal_del = c(CG_female_paternal_SNV_counts$all, CG_female_paternal_SNV_counts$syn, CG_female_paternal_SNV_counts$mis, CG_female_paternal_SNV_counts$lof),
                              father_of_female = c(CG_father_of_female_SNV_counts$all, CG_father_of_female_SNV_counts$syn, CG_father_of_female_SNV_counts$mis, CG_father_of_female_SNV_counts$lof),
                              male_proband_maternal_del = c(CG_male_maternal_SNV_counts$all, CG_male_maternal_SNV_counts$syn, CG_male_maternal_SNV_counts$mis, CG_male_maternal_SNV_counts$lof),
                              mother_of_male = c(CG_mother_of_male_SNV_counts$all, CG_mother_of_male_SNV_counts$syn, CG_mother_of_male_SNV_counts$mis, CG_mother_of_male_SNV_counts$lof),
                              female_proband_maternal_del = c(CG_female_maternal_SNV_counts$all, CG_female_maternal_SNV_counts$syn, CG_female_maternal_SNV_counts$mis, CG_female_maternal_SNV_counts$lof),
                              mother_of_female = c(CG_mother_of_female_SNV_counts$all, CG_mother_of_female_SNV_counts$syn, CG_mother_of_female_SNV_counts$mis, CG_mother_of_female_SNV_counts$lof))

CG_fullSNV_comparison_df <- as.data.frame(CG_fullSNV_comparison)

CG_full_indiv_count_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                       male_proband_paternal_del = c(CG_male_paternal_counts$all, CG_male_paternal_counts$syn, CG_male_paternal_counts$mis, CG_male_paternal_counts$lof),
                                       father_of_male = c(CG_father_of_male_counts$all, CG_father_of_male_counts$syn, CG_father_of_male_counts$mis, CG_father_of_male_counts$lof),
                                       female_proband_paternal_del = c(CG_female_paternal_counts$all, CG_female_paternal_counts$syn, CG_female_paternal_counts$mis, CG_female_paternal_counts$lof),
                                       father_of_female = c(CG_father_of_female_counts$all, CG_father_of_female_counts$syn, CG_father_of_female_counts$mis, CG_father_of_female_counts$lof),
                                       male_proband_maternal_del = c(CG_male_maternal_counts$all, CG_male_maternal_counts$syn, CG_male_maternal_counts$mis, CG_male_maternal_counts$lof),
                                       mother_of_male = c(CG_mother_of_male_counts$all, CG_mother_of_male_counts$syn, CG_mother_of_male_counts$mis, CG_mother_of_male_counts$lof),
                                       female_proband_maternal_del = c(CG_female_maternal_counts$all, CG_female_maternal_counts$syn, CG_female_maternal_counts$mis, CG_female_maternal_counts$lof),
                                       mother_of_female = c(CG_mother_of_female_counts$all, CG_mother_of_female_counts$syn, CG_mother_of_female_counts$mis, CG_mother_of_female_counts$lof))

CG_full_indiv_count_comparison_df <- as.data.frame(CG_full_indiv_count_comparison)

MSSNG_fullSNV_comparison_df <- cbind("Variant Type"=ILMN_fullSNV_comparison_df[, 1], ILMN_fullSNV_comparison_df[, -1] + CG_fullSNV_comparison_df[, -1])
MSSNG_full_indiv_count_comparison_df <- cbind("Variant Type"=ILMN_full_indiv_count_comparison_df[, 1], ILMN_full_indiv_count_comparison_df[, -1] + CG_full_indiv_count_comparison_df[, -1])

# With only probands grouped
MSSNG_SNV_comparison_df <- data.frame("Variant Type"=MSSNG_fullSNV_comparison_df[, 1], 
                                 "proband_paternal_del"=MSSNG_fullSNV_comparison_df[, 2]+MSSNG_fullSNV_comparison_df[, 4],
                                 "father"=MSSNG_fullSNV_comparison_df[, 3]+MSSNG_fullSNV_comparison_df[, 5],
                                 "proband_maternal_del"=MSSNG_fullSNV_comparison_df[, 6]+MSSNG_fullSNV_comparison_df[, 8],
                                 "mother"=MSSNG_fullSNV_comparison_df[, 7]+MSSNG_fullSNV_comparison_df[, 9])

MSSNG_indiv_count_comparison_df <- data.frame("Variant Type"=MSSNG_full_indiv_count_comparison_df[, 1], 
                                         "proband_paternal_del"=MSSNG_full_indiv_count_comparison_df[, 2]+MSSNG_full_indiv_count_comparison_df[, 4],
                                         "father"=MSSNG_full_indiv_count_comparison_df[, 3]+MSSNG_full_indiv_count_comparison_df[, 5],
                                         "proband_maternal_del"=MSSNG_full_indiv_count_comparison_df[, 6]+MSSNG_full_indiv_count_comparison_df[, 8],
                                         "mother"=MSSNG_full_indiv_count_comparison_df[, 7]+MSSNG_full_indiv_count_comparison_df[, 9])

# Full grouped (all probands vs all parents)
MSSNG_groupedSNV_comparison_df <- data.frame("Variant Type"=MSSNG_SNV_comparison_df[, 1], "Probands"=MSSNG_SNV_comparison_df[, 2] + MSSNG_SNV_comparison_df[, 4], "Parents"=MSSNG_SNV_comparison_df[, 3] + MSSNG_SNV_comparison_df[, 5])
MSSNG_grouped_indiv_count_comparison_df <- data.frame("Variant Type"=MSSNG_indiv_count_comparison_df[, 1], "Probands"=MSSNG_indiv_count_comparison_df[, 2] + MSSNG_indiv_count_comparison_df[, 4], "Parents"=MSSNG_indiv_count_comparison_df[, 3] + MSSNG_indiv_count_comparison_df[, 5])

##### SSC Data Counts ######
SSC_sex_IDs <- Get_IDs_By_Sex("SSC_metadata.tsv")
SSC_male_IDs <- SSC_sex_IDs[[1]]
SSC_father_of_male_IDs <- SSC_sex_IDs[[2]]
SSC_mother_of_male_IDs <- SSC_sex_IDs[[3]]

SSC_female_IDs <- SSC_sex_IDs[[4]]
SSC_father_of_female_IDs <- SSC_sex_IDs[[5]]
SSC_mother_of_female_IDs <- SSC_sex_IDs[[6]]

SSC_male_probandSNVs_in_paternal_CNVs <- SSC_probandSNVs_in_paternalCNVs[which(SSC_probandSNVs_in_paternalCNVs$X.Sample %in% SSC_male_IDs), ]
SSC_male_probandSNVs_in_maternal_CNVs <- SSC_probandSNVs_in_maternalCNVs[which(SSC_probandSNVs_in_maternalCNVs$X.Sample %in% SSC_male_IDs), ]
SSC_father_of_male_SNVs_in_paternal_CNVs <- SSC_fatherSNVs_in_paternalCNVs[which(SSC_fatherSNVs_in_paternalCNVs$X.Sample %in% SSC_father_of_male_IDs), ]
SSC_mother_of_male_SNVs_in_maternal_CNVs <- SSC_motherSNVs_in_maternalCNVs[which(SSC_motherSNVs_in_maternalCNVs$X.Sample %in% SSC_mother_of_male_IDs), ]

SSC_female_probandSNVs_in_paternal_CNVs <- SSC_probandSNVs_in_paternalCNVs[which(SSC_probandSNVs_in_paternalCNVs$X.Sample %in% SSC_female_IDs), ]
SSC_female_probandSNVs_in_maternal_CNVs <- SSC_probandSNVs_in_maternalCNVs[which(SSC_probandSNVs_in_maternalCNVs$X.Sample %in% SSC_female_IDs), ]
SSC_father_of_female_SNVs_in_paternal_CNVs <- SSC_fatherSNVs_in_paternalCNVs[which(SSC_fatherSNVs_in_paternalCNVs$X.Sample %in% SSC_father_of_female_IDs), ]
SSC_mother_of_female_SNVs_in_maternal_CNVs <- SSC_motherSNVs_in_maternalCNVs[which(SSC_motherSNVs_in_maternalCNVs$X.Sample %in% SSC_mother_of_female_IDs), ]

SSC_male_paternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_male_probandSNVs_in_paternal_CNVs)
SSC_male_paternal_SNV_counts <- SSC_male_paternal_both_counts[[1]]
SSC_male_paternal_counts <- SSC_male_paternal_both_counts[[2]]

SSC_male_maternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_male_probandSNVs_in_maternal_CNVs)
SSC_male_maternal_SNV_counts <- SSC_male_maternal_both_counts[[1]]
SSC_male_maternal_counts <- SSC_male_maternal_both_counts[[2]]

SSC_father_of_male_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_father_of_male_SNVs_in_paternal_CNVs)
SSC_father_of_male_SNV_counts <- SSC_father_of_male_both_counts[[1]]
SSC_father_of_male_counts <- SSC_father_of_male_both_counts[[2]]

SSC_mother_of_male_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_mother_of_male_SNVs_in_maternal_CNVs)
SSC_mother_of_male_SNV_counts <- SSC_mother_of_male_both_counts[[1]]
SSC_mother_of_male_counts <- SSC_mother_of_male_both_counts[[2]]

SSC_female_paternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_female_probandSNVs_in_paternal_CNVs)
SSC_female_paternal_SNV_counts <- SSC_female_paternal_both_counts[[1]]
SSC_female_paternal_counts <- SSC_female_paternal_both_counts[[2]]

SSC_female_maternal_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_female_probandSNVs_in_maternal_CNVs)
SSC_female_maternal_SNV_counts <- SSC_female_maternal_both_counts[[1]]
SSC_female_maternal_counts <- SSC_female_maternal_both_counts[[2]]

SSC_father_of_female_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_father_of_female_SNVs_in_paternal_CNVs)
SSC_father_of_female_SNV_counts <- SSC_father_of_female_both_counts[[1]]
SSC_father_of_female_counts <- SSC_father_of_female_both_counts[[2]]

SSC_mother_of_female_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_mother_of_female_SNVs_in_maternal_CNVs)
SSC_mother_of_female_SNV_counts <- SSC_mother_of_female_both_counts[[1]]
SSC_mother_of_female_counts <- SSC_mother_of_female_both_counts[[2]]

SSC_fullSNV_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                               male_proband_paternal_del = c(SSC_male_paternal_SNV_counts$all, SSC_male_paternal_SNV_counts$syn, SSC_male_paternal_SNV_counts$mis, SSC_male_paternal_SNV_counts$lof),
                               father_of_male = c(SSC_father_of_male_SNV_counts$all, SSC_father_of_male_SNV_counts$syn, SSC_father_of_male_SNV_counts$mis, SSC_father_of_male_SNV_counts$lof),
                               female_proband_paternal_del = c(SSC_female_paternal_SNV_counts$all, SSC_female_paternal_SNV_counts$syn, SSC_female_paternal_SNV_counts$mis, SSC_female_paternal_SNV_counts$lof),
                               father_of_female = c(SSC_father_of_female_SNV_counts$all, SSC_father_of_female_SNV_counts$syn, SSC_father_of_female_SNV_counts$mis, SSC_father_of_female_SNV_counts$lof),
                               male_proband_maternal_del = c(SSC_male_maternal_SNV_counts$all, SSC_male_maternal_SNV_counts$syn, SSC_male_maternal_SNV_counts$mis, SSC_male_maternal_SNV_counts$lof),
                               mother_of_male = c(SSC_mother_of_male_SNV_counts$all, SSC_mother_of_male_SNV_counts$syn, SSC_mother_of_male_SNV_counts$mis, SSC_mother_of_male_SNV_counts$lof),
                               female_proband_maternal_del = c(SSC_female_maternal_SNV_counts$all, SSC_female_maternal_SNV_counts$syn, SSC_female_maternal_SNV_counts$mis, SSC_female_maternal_SNV_counts$lof),
                               mother_of_female = c(SSC_mother_of_female_SNV_counts$all, SSC_mother_of_female_SNV_counts$syn, SSC_mother_of_female_SNV_counts$mis, SSC_mother_of_female_SNV_counts$lof))

SSC_fullSNV_comparison_df <- as.data.frame(SSC_fullSNV_comparison)

SSC_full_indiv_count_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                        male_proband_paternal_del = c(SSC_male_paternal_counts$all, SSC_male_paternal_counts$syn, SSC_male_paternal_counts$mis, SSC_male_paternal_counts$lof),
                                        father_of_male = c(SSC_father_of_male_counts$all, SSC_father_of_male_counts$syn, SSC_father_of_male_counts$mis, SSC_father_of_male_counts$lof),
                                        female_proband_paternal_del = c(SSC_female_paternal_counts$all, SSC_female_paternal_counts$syn, SSC_female_paternal_counts$mis, SSC_female_paternal_counts$lof),
                                        father_of_female = c(SSC_father_of_female_counts$all, SSC_father_of_female_counts$syn, SSC_father_of_female_counts$mis, SSC_father_of_female_counts$lof),
                                        male_proband_maternal_del = c(SSC_male_maternal_counts$all, SSC_male_maternal_counts$syn, SSC_male_maternal_counts$mis, SSC_male_maternal_counts$lof),
                                        mother_of_male = c(SSC_mother_of_male_counts$all, SSC_mother_of_male_counts$syn, SSC_mother_of_male_counts$mis, SSC_mother_of_male_counts$lof),
                                        female_proband_maternal_del = c(SSC_female_maternal_counts$all, SSC_female_maternal_counts$syn, SSC_female_maternal_counts$mis, SSC_female_maternal_counts$lof),
                                        mother_of_female = c(SSC_mother_of_female_counts$all, SSC_mother_of_female_counts$syn, SSC_mother_of_female_counts$mis, SSC_mother_of_female_counts$lof))

SSC_full_indiv_count_comparison_df <- as.data.frame(SSC_full_indiv_count_comparison)

# With only probands grouped
SSC_SNV_comparison_df <- data.frame("Variant Type"=SSC_fullSNV_comparison_df[, 1], 
                               "proband_paternal_del"=SSC_fullSNV_comparison_df[, 2]+SSC_fullSNV_comparison_df[, 4],
                               "father"=SSC_fullSNV_comparison_df[, 3]+SSC_fullSNV_comparison_df[, 5],
                               "proband_maternal_del"=SSC_fullSNV_comparison_df[, 6]+SSC_fullSNV_comparison_df[, 8],
                               "mother"=SSC_fullSNV_comparison_df[, 7]+SSC_fullSNV_comparison_df[, 9])

SSC_indiv_count_comparison_df <- data.frame("Variant Type"=SSC_full_indiv_count_comparison_df[, 1], 
                                       "proband_paternal_del"=SSC_full_indiv_count_comparison_df[, 2]+SSC_full_indiv_count_comparison_df[, 4],
                                       "father"=SSC_full_indiv_count_comparison_df[, 3]+SSC_full_indiv_count_comparison_df[, 5],
                                       "proband_maternal_del"=SSC_full_indiv_count_comparison_df[, 6]+SSC_full_indiv_count_comparison_df[, 8],
                                       "mother"=SSC_full_indiv_count_comparison_df[, 7]+SSC_full_indiv_count_comparison_df[, 9])

# Full grouped (all probands vs all parents)
SSC_groupedSNV_comparison_df <- data.frame("Variant Type"=SSC_SNV_comparison_df[, 1], "Probands"=SSC_SNV_comparison_df[, 2] + SSC_SNV_comparison_df[, 4], "Parents"=SSC_SNV_comparison_df[, 3] + SSC_SNV_comparison_df[, 5])
SSC_grouped_indiv_count_comparison_df <- data.frame("Variant Type"=SSC_indiv_count_comparison_df[, 1], "Probands"=SSC_indiv_count_comparison_df[, 2] + SSC_indiv_count_comparison_df[, 4], "Parents"=SSC_indiv_count_comparison_df[, 3] + SSC_indiv_count_comparison_df[, 5])


##### SSC DATA - Probands vs. Unaffected Siblings ######

SSC_unaffectedSibling_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_unaffectedSiblingsSNVs)
SSC_unaffectedSibling_SNV_counts <- SSC_unaffectedSibling_both_counts[[1]]
SSC_unaffectedSibling_counts <- SSC_unaffectedSibling_both_counts[[2]]

SSC_probandUN_both_counts <- Get_SNV_Count_By_Var_Type(loaded_data = SSC_probandUNSNVs)
SSC_probandUN_SNV_counts <- SSC_probandUN_both_counts[[1]]
SSC_probandUN_counts <- SSC_probandUN_both_counts[[2]]


SSC_SNV_sibling_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                                         proband_del = c(SSC_probandUN_SNV_counts$all, SSC_probandUN_SNV_counts$syn, SSC_probandUN_SNV_counts$mis, SSC_probandUN_SNV_counts$lof),
                                                         unaffectedSibling_del = c(SSC_unaffectedSibling_SNV_counts$all, SSC_unaffectedSibling_SNV_counts$syn, SSC_unaffectedSibling_SNV_counts$mis, SSC_unaffectedSibling_SNV_counts$lof))

SSC_SNV_sibling_comparison_df <- as.data.frame(SSC_SNV_sibling_comparison)

SSC_SNV_sibling_comparison_df_norm <- SSC_SNV_sibling_comparison_df
SSC_SNV_sibling_comparison_df_norm[, 2] <- SSC_SNV_sibling_comparison_df_norm[, 2]/(SSC_probandUN_exonic_size/1000000)
SSC_SNV_sibling_comparison_df_norm[, 3] <- SSC_SNV_sibling_comparison_df_norm[, 3]/(SSC_unaffectedSiblings_exonic_size/1000000)


SSC_SNV_sibling_indiv_count_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                   proband_del = c(SSC_probandUN_counts$all, SSC_probandUN_counts$syn, SSC_probandUN_counts$mis, SSC_probandUN_counts$lof),
                                   unaffectedSibling_del = c(SSC_unaffectedSibling_counts$all, SSC_unaffectedSibling_counts$syn, SSC_unaffectedSibling_counts$mis, SSC_unaffectedSibling_counts$lof))

SSC_SNV_sibling_indiv_count_comparison_df <- as.data.frame(SSC_SNV_sibling_indiv_count_comparison)

