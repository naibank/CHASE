#Required Packages
#BiocManager::install("GenomicRanges")
#BiocManager::install("Repitools")

library(GenomicRanges)
library(Repitools)
library(yaml)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE/Faraz")

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

Get_SNVs_Parent_Comparison <- function(file_path, CNV_freq_filter=1, SNV_freq_filter=1) {
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
          childSNVs.in.maternalCNVs <- rbind(childSNVs.in.maternalCNVs, childSNVs[unique(olap@from), ])
        }
        
        if(nrow(motherCNVs) > 0 & nrow(motherSNVs) > 0){
          names(motherSNVs)[11] <- "SampleData"
          motherCNVs <- GRanges(motherCNVs$chrAnn,
                                IRanges(motherCNVs$STARTAnn,
                                        motherCNVs$ENDAnn), "*")
          motherSNVs <- motherSNVs[,names(motherSNVs)[names(motherSNVs) != "Inheritance.SNV"]]
          motherSNVs.g <- GRanges(motherSNVs$CHROM, IRanges(motherSNVs$start, motherSNVs$end), "*")
          olap <- findOverlaps(motherSNVs.g, motherCNVs)
          motherSNVs.in.maternalCNVs <- rbind(motherSNVs.in.maternalCNVs, motherSNVs[unique(olap@from), ])
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
          childSNVs.in.paternalCNVs <- rbind(childSNVs.in.paternalCNVs, childSNVs[unique(olap@from), ])
        }
        
        if(nrow(fatherCNVs) > 0 & nrow(fatherSNVs) > 0){
          names(fatherSNVs)[11] <- "SampleData"
          fatherCNVs <- GRanges(fatherCNVs$chrAnn,
                                IRanges(fatherCNVs$STARTAnn,
                                        fatherCNVs$ENDAnn), "*")
          fatherSNVs <- fatherSNVs[,names(fatherSNVs)[names(fatherSNVs) != "Inheritance.SNV"]]
          
          fatherSNVs.g <- GRanges(fatherSNVs$CHROM, IRanges(fatherSNVs$start, fatherSNVs$end), "*")
          olap <- findOverlaps(fatherSNVs.g, fatherCNVs)
          fatherSNVs.in.paternalCNVs <- rbind(fatherSNVs.in.paternalCNVs, fatherSNVs[unique(olap@from), ])
        }
      }
    }
    message(i)
  }
  
  return (list(childSNVs.in.paternalCNVs, 
               childSNVs.in.maternalCNVs, 
               fatherSNVs.in.paternalCNVs, 
               motherSNVs.in.maternalCNVs))
}

Get_SNVs_Sibling_Comparison <- function(file_path, CNV_freq_filter=1, SNV_freq_filter=1) {
  data <- yaml::yaml.load_file(file_path)
  
  childSNVs_combined <- data.frame()
  SNVs_to_add <- data.frame()
  
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
      
      if(nrow(childSNVs) > 0){
        names(childSNVs)[11] <- "SampleData"
        
        childSNVs_g <- GRanges(childSNVs$CHROM, IRanges(childSNVs$start, childSNVs$end), "*")
        olap <- findOverlaps(childSNVs_g, childCNVs_g)
        SNVs_to_add <- childSNVs[unique(olap@from), ]
        SNVs_to_add$cnvInheritance <- childCNVs$Inheritance[olap@to]
        childSNVs_combined <- rbind(childSNVs_combined, SNVs_to_add)
      }
    }
    message(i)
  }
  
  # Exclude SNVs inherited from the deletion transmitting parent
  childSNVs_combined <- subset(childSNVs_combined, 
                               childSNVs_combined$inheritance != tolower(childSNVs_combined$cnvInheritance))
  
  
  return (childSNVs_combined)
}

##### ILMN DATA ##### (1111 CH events)

ILMN_SNVs_data <- Get_SNVs_Parent_Comparison("MSSNG_ILMN_CH_Data10P_Bank.yaml", CNV_freq_filter=0.01)

ILMN_probandSNVs_in_paternalCNVs <- ILMN_SNVs_data[[1]]
ILMN_probandSNVs_in_maternalCNVs <- ILMN_SNVs_data[[2]]
ILMN_fatherSNVs_in_maternalCNVs <- ILMN_SNVs_data[[3]]
ILMN_motherSNVs_in_paternalCNVs <- ILMN_SNVs_data[[4]]

# Write SNVs to tsv
write.table(ILMN_probandSNVs_in_paternalCNVs, "./SNV Data/ILMN_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(ILMN_probandSNVs_in_maternalCNVs, "./SNV Data/ILMN_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(ILMN_fatherSNVs_in_maternalCNVs, "./SNV Data/ILMN_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(ILMN_motherSNVs_in_paternalCNVs, "./SNV Data/ILMN_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

##### CG DATA ##### (120 CH events)

CG_SNVs_data <- Get_SNVs_Parent_Comparison("MSSNG_CG_CH_Data10P_Bank.yaml", CNV_freq_filter=0.01)

CG_probandSNVs_in_paternalCNVs <- CG_SNVs_data[[1]]
CG_probandSNVs_in_maternalCNVs <- CG_SNVs_data[[2]]
CG_fatherSNVs_in_maternalCNVs <- CG_SNVs_data[[3]]
CG_motherSNVs_in_paternalCNVs <- CG_SNVs_data[[4]]

# Write SNVs to tsv
write.table(CG_probandSNVs_in_paternalCNVs, "./SNV Data/CG_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(CG_probandSNVs_in_maternalCNVs, "./SNV Data/CG_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(CG_fatherSNVs_in_maternalCNVs, "./SNV Data/CG_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(CG_motherSNVs_in_paternalCNVs, "./SNV Data/CG_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

##### SSC DATA ##### (1510 CH events)

SSC_SNVs_data <- Get_SNVs_Parent_Comparison("SSC_CH_Data10P_Bank.yaml", CNV_freq_filter=0.01)

SSC_probandSNVs_in_paternalCNVs <- SSC_SNVs_data[[1]]
SSC_probandSNVs_in_maternalCNVs <- SSC_SNVs_data[[2]]
SSC_fatherSNVs_in_maternalCNVs <- SSC_SNVs_data[[3]]
SSC_motherSNVs_in_paternalCNVs <- SSC_SNVs_data[[4]]

# Write SNVs to tsv
write.table(SSC_probandSNVs_in_paternalCNVs, "./SNV Data/SSC_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_probandSNVs_in_maternalCNVs, "./SNV Data/SSC_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_fatherSNVs_in_maternalCNVs, "./SNV Data/SSC_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_motherSNVs_in_paternalCNVs, "./SNV Data/SSC_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

##### SSC DATA - Proband vs. Unaffected Siblings ##### (X CH events)

SSC_probandUNSNVs <- Get_SNVs_Sibling_Comparison("SSC_CH_Data10P_Bank.yaml", CNV_freq_filter=0.01)
SSC_unaffectedSiblingsSNVs <- Get_SNVs_Sibling_Comparison("SSC_CH.unaffectedSiblings_Data10P_Bank.yaml", CNV_freq_filter=0.01)

# Write SNVs to tsv
write.table(SSC_probandUNSNVs, "./SNV Data/SSC_ProbandUNSNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_unaffectedSiblingsSNVs, "./SNV Data/SSC_UnaffectedSiblingsSNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

###########################################################################################################################################################

# SNV COUNT BY VARIANT TYPE

Get_SNV_Count_By_Var_Type <- function(file_path) {
  # Read SNVs tsv
  SNVs_data <- read.delim(file_path, stringsAsFactors = F)
  
  # All variants
  SNVs_all <- length(SNVs_data$X.Sample[!duplicated(paste(SNVs_data$X.Sample, SNVs_data$X.id))])                
  SNvs_all_count <- length(unique(SNVs_data$X.Sample))
  
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
  
  return (list(all=SNVs_all, syn=SNVs_syn, mis=SNVs_mis, lof=SNVs_lof))
}

##### ILMN DATA ######

ILMN_proband_paternal_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/ILMN_ProbandSNVs_in_PaternalCNVs.tsv")
ILMN_proband_maternal_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/ILMN_ProbandSNVs_in_MaternalCNVs.tsv")
ILMN_father_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/ILMN_FatherSNVs_in_PaternalCNVs.tsv")
ILMN_mother_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/ILMN_MotherSNVs_in_MaternalCNVs.tsv")

ILMN_SNV_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                      proband_paternal_del = c(ILMN_proband_paternal_del_SNVs$all, ILMN_proband_paternal_del_SNVs$syn, ILMN_proband_paternal_del_SNVs$mis, ILMN_proband_paternal_del_SNVs$lof),
                      proband_maternal_del = c(ILMN_proband_maternal_del_SNVs$all, ILMN_proband_maternal_del_SNVs$syn, ILMN_proband_maternal_del_SNVs$mis, ILMN_proband_maternal_del_SNVs$lof), 
                      father = c(ILMN_father_del_SNVs$all, ILMN_father_del_SNVs$syn, ILMN_father_del_SNVs$mis, ILMN_father_del_SNVs$lof),
                      mother = c(ILMN_mother_del_SNVs$all, ILMN_mother_del_SNVs$syn, ILMN_mother_del_SNVs$mis, ILMN_mother_del_SNVs$lof))

ILMN_SNV_comparison_df <- as.data.frame(ILMN_SNV_comparison)

##### CG DATA ######

CG_proband_paternal_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/CG_ProbandSNVs_in_PaternalCNVs.tsv")
CG_proband_maternal_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/CG_ProbandSNVs_in_MaternalCNVs.tsv")
CG_father_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/CG_FatherSNVs_in_PaternalCNVs.tsv")
CG_mother_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/CG_MotherSNVs_in_MaternalCNVs.tsv")

CG_SNV_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                            proband_paternal_del = c(CG_proband_paternal_del_SNVs$all, CG_proband_paternal_del_SNVs$syn, CG_proband_paternal_del_SNVs$mis, CG_proband_paternal_del_SNVs$lof),
                            proband_maternal_del = c(CG_proband_maternal_del_SNVs$all, CG_proband_maternal_del_SNVs$syn, CG_proband_maternal_del_SNVs$mis, CG_proband_maternal_del_SNVs$lof), 
                            father = c(CG_father_del_SNVs$all, CG_father_del_SNVs$syn, CG_father_del_SNVs$mis, CG_father_del_SNVs$lof),
                            mother = c(CG_mother_del_SNVs$all, CG_mother_del_SNVs$syn, CG_mother_del_SNVs$mis, CG_mother_del_SNVs$lof))

CG_SNV_comparison_df <- as.data.frame(CG_SNV_comparison)

##### MSSNG COMBINED #####
MSSNG_SNV_comparison_df <- cbind("Variant Type"=ILMN_SNV_comparison_df[, 1], ILMN_SNV_comparison_df[, -1] + CG_SNV_comparison_df[, -1])

##### SSC DATA ######

SSC_proband_paternal_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/SSC_ProbandSNVs_in_PaternalCNVs.tsv")
SSC_proband_maternal_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/SSC_ProbandSNVs_in_MaternalCNVs.tsv")
SSC_father_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/SSC_FatherSNVs_in_PaternalCNVs.tsv")
SSC_mother_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/SSC_MotherSNVs_in_MaternalCNVs.tsv")

SSC_SNV_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                          proband_paternal_del = c(SSC_proband_paternal_del_SNVs$all, SSC_proband_paternal_del_SNVs$syn, SSC_proband_paternal_del_SNVs$mis, SSC_proband_paternal_del_SNVs$lof),
                          proband_maternal_del = c(SSC_proband_maternal_del_SNVs$all, SSC_proband_maternal_del_SNVs$syn, SSC_proband_maternal_del_SNVs$mis, SSC_proband_maternal_del_SNVs$lof), 
                          father = c(SSC_father_del_SNVs$all, SSC_father_del_SNVs$syn, SSC_father_del_SNVs$mis, SSC_father_del_SNVs$lof),
                          mother = c(SSC_mother_del_SNVs$all, SSC_mother_del_SNVs$syn, SSC_mother_del_SNVs$mis, SSC_mother_del_SNVs$lof))

SSC_SNV_comparison_df <- as.data.frame(SSC_SNV_comparison)

SSC_unaffectedSibling_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/SSC_UnaffectedSiblingsSNVs.tsv")
SSC_probandUN_del_SNVs <- Get_SNV_Count_By_Var_Type("./SNV Data/SSC_ProbandUNSNVs.tsv")

SSC_SNV_sibling_comparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                                                         proband_del = c(SSC_probandUN_del_SNVs$all, SSC_probandUN_del_SNVs$syn, SSC_probandUN_del_SNVs$mis, SSC_probandUN_del_SNVs$lof),
                                                         unaffectedSibling_del = c(SSC_unaffectedSibling_del_SNVs$all, SSC_unaffectedSibling_del_SNVs$syn, SSC_unaffectedSibling_del_SNVs$mis, SSC_unaffectedSibling_del_SNVs$lof))

SSC_SNV_sibling_comparison_df <- as.data.frame(SSC_SNV_sibling_comparison)

