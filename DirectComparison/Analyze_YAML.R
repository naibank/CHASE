#Required Packages
#BiocManager::install("GenomicRanges")
#BiocManager::install("Repitools")

library(GenomicRanges)
library(Repitools)
library(yaml)

setwd("C:/Users/acale/Desktop/Fall Co-op 2021/CHASE Files/MSSNG")

##### ILMN DATA ##### (1111 CH events)

ILMNdata <- yaml::yaml.load_file("MSSNG_ILMN_CH_Data10P_Bank.yaml")

probandSNVs.in.maternalCNVs <- data.frame()
probandSNVs.in.paternalCNVs <- data.frame()
motherSNVs.in.maternalCNVs <- data.frame()
fatherSNVs.in.paternalCNVs <- data.frame()

maternalCount = 0
paternalCount = 0

motherCount = 0
fatherCount = 0
probandBothPaternalMaternalCount = 0

for(i in 1:length(ILMNdata)){
  probandCNVs <- data.frame(ILMNdata[[i]]$CNVs)
  
  motherCount = motherCount + as.numeric("Maternal" %in% probandCNVs$Inheritance)
  fatherCount = fatherCount + as.numeric("Paternal" %in% probandCNVs$Inheritance)
  probandBothPaternalMaternalCount = probandBothPaternalMaternalCount + as.numeric("Paternal" %in% probandCNVs$Inheritance &
                                                                                     "Maternal" %in% probandCNVs$Inheritance)
  
  probandSNVs <- data.frame(ILMNdata[[i]]$SNVs)
  fatherCNVs <- data.frame(ILMNdata[[i]]$Father.CNV)
  fatherSNVs <- data.frame(ILMNdata[[i]]$Father.SNV)
  motherCNVs <- data.frame(ILMNdata[[i]]$Mother.CNV)
  motherSNVs <- data.frame(ILMNdata[[i]]$Mother.SNV)
  
  probandSNVs <- probandSNVs[which(probandSNVs$freq_max < 0.1), ]
  fatherSNVs <- fatherSNVs[which(fatherSNVs$freq_max < 0.1), ]
  motherSNVs <- motherSNVs[which(motherSNVs$freq_max < 0.1), ]
  
  if(!("Paternal" %in% probandCNVs$Inheritance & "Maternal" %in% probandCNVs$Inheritance)){
    if(sum(probandCNVs$Inheritance == "Maternal") > 0){ # & nrow(motherCNVs) > 0){
      maternalCount = maternalCount + 1 
      
      message(sprintf("%s %s %s", i, motherCount-probandBothPaternalMaternalCount, maternalCount))
      #motherCNVs <- GRanges(motherCNVs$chrAnn, IRanges(motherCNVs$STARTAnn, motherCNVs$ENDAnn), "*")
      proband.motherCNVs <- GRanges(probandCNVs$chrAnn[probandCNVs$Inheritance == "Maternal"],
                                   IRanges(probandCNVs$STARTAnn[probandCNVs$Inheritance == "Maternal"],
                                           probandCNVs$ENDAnn[probandCNVs$Inheritance == "Maternal"]), "*")
      motherCNVs <- GRanges(motherCNVs$chrAnn,
                            IRanges(motherCNVs$STARTAnn,
                                    motherCNVs$ENDAnn), "*")
    
      if(nrow(probandSNVs) > 0){
        names(probandSNVs)[11] <- "SampleData"
        
        probandSNVs.g <- GRanges(probandSNVs$CHROM, IRanges(probandSNVs$start, probandSNVs$end), "*")
        olap <- findOverlaps(probandSNVs.g, proband.motherCNVs)
        probandSNVs.in.maternalCNVs <- rbind(probandSNVs.in.maternalCNVs, probandSNVs[unique(olap@from), ])
      }
    
      if(nrow(motherSNVs) > 0){
        names(motherSNVs)[11] <- "SampleData"
        motherSNVs <- motherSNVs[,names(motherSNVs)[names(motherSNVs) != "Inheritance.SNV"]]
        motherSNVs.g <- GRanges(motherSNVs$CHROM, IRanges(motherSNVs$start, motherSNVs$end), "*")
        olap <- findOverlaps(motherSNVs.g, motherCNVs)
        motherSNVs.in.maternalCNVs <- rbind(motherSNVs.in.maternalCNVs, motherSNVs[unique(olap@from), ])
      }
    }
  
    if(sum(probandCNVs$Inheritance == "Paternal") > 0){ #} & nrow(fatherCNVs) > 0){
      paternalCount = paternalCount + 1
      # fatherCNVs <- GRanges(fatherCNVs$chrAnn, IRanges(fatherCNVs$STARTAnn, fatherCNVs$ENDAnn), "*")
      proband.fatherCNVs <- GRanges(probandCNVs$chrAnn[probandCNVs$Inheritance == "Paternal"],
                                   IRanges(probandCNVs$STARTAnn[probandCNVs$Inheritance == "Paternal"],
                                           probandCNVs$ENDAnn[probandCNVs$Inheritance == "Paternal"]), "*")
      fatherCNVs <- GRanges(fatherCNVs$chrAnn,
                            IRanges(fatherCNVs$STARTAnn,
                                    fatherCNVs$ENDAnn), "*")
      if(nrow(probandSNVs) > 0){
        names(probandSNVs)[11] <- "SampleData"
      
        probandSNVs.g <- GRanges(probandSNVs$CHROM, IRanges(probandSNVs$start, probandSNVs$end), "*")
        olap <- findOverlaps(probandSNVs.g, proband.fatherCNVs)
        probandSNVs.in.paternalCNVs <- rbind(probandSNVs.in.paternalCNVs, probandSNVs[unique(olap@from), ])
      }
    
      if(nrow(fatherSNVs) > 0){
        names(fatherSNVs)[11] <- "SampleData"
        fatherSNVs <- fatherSNVs[,names(fatherSNVs)[names(fatherSNVs) != "Inheritance.SNV"]]
      
        fatherSNVs.g <- GRanges(fatherSNVs$CHROM, IRanges(fatherSNVs$start, fatherSNVs$end), "*")
        olap <- findOverlaps(fatherSNVs.g, fatherCNVs)
        fatherSNVs.in.paternalCNVs <- rbind(fatherSNVs.in.paternalCNVs, fatherSNVs[unique(olap@from), ])
      }
    }
  }
  message(i)
}

#Remove probands that inherited CNVs from both mother and father
# probandSNVs.in.paternalCNVsNA <- probandSNVs.in.paternalCNVs[ !is.na(probandSNVs.in.paternalCNVs$Sample.ID), ]
# probandSNVs.in.maternalCNVsNA <- probandSNVs.in.maternalCNVs[ !is.na(probandSNVs.in.maternalCNVs$Sample.ID), ]
# 
# i = 0
# 
# probandCH_SNVs.duplicate <- vector("character", length = 1)
# 
# for(j in 1:NROW(probandSNVs.in.paternalCNVsNA)){
#   for(k in 1:NROW(probandSNVs.in.maternalCNVsNA)){
#     if(probandSNVs.in.paternalCNVsNA$Sample.ID[j] == probandSNVs.in.maternalCNVsNA$Sample.ID[k]){
#       i = i + 1
#       probandCH_SNVs.duplicate[i] <- probandSNVs.in.paternalCNVsNA$Sample.ID[j]
#     }
#   }
# }
# 
# probandSNVs.in.paternalCNVs <- probandSNVs.in.paternalCNVs[ !probandSNVs.in.paternalCNVs$Sample.ID %in% probandCH_SNVs.duplicate, ]
# probandSNVs.in.maternalCNVs <- probandSNVs.in.maternalCNVs[ !probandSNVs.in.maternalCNVs$Sample.ID %in% probandCH_SNVs.duplicate, ]

#Write SNVs to tsv
write.table(probandSNVs.in.paternalCNVs, "ILMN_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(probandSNVs.in.maternalCNVs, "ILMN_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(fatherSNVs.in.paternalCNVs, "ILMN_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(motherSNVs.in.maternalCNVs, "ILMN_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)


##### CG DATA ##### (120 CH events)

CGdata <- yaml::yaml.load_file("MSSNG_CG_CH_Data10P_Bank.yaml")

probandSNVs.CG.in.maternalCNVs <- data.frame()
probandSNVs.CG.in.paternalCNVs <- data.frame()
motherSNVs.CG.in.maternalCNVs <- data.frame()
fatherSNVs.CG.in.paternalCNVs <- data.frame()

maternalCount.CG = 0
paternalCount.CG = 0

motherCount.CG = 0
fatherCount.CG = 0
probandBothPaternalMaternalCount.CG = 0

for(i in 1:length(CGdata)){
  probandCNVs.CG <- data.frame(CGdata[[i]]$CNVs)
  
  motherCount.CG = motherCount.CG + as.numeric("Maternal" %in% probandCNVs.CG$Inheritance)
  fatherCount.CG = fatherCount.CG + as.numeric("Paternal" %in% probandCNVs.CG$Inheritance)
  probandBothPaternalMaternalCount.CG = probandBothPaternalMaternalCount.CG + as.numeric("Paternal" %in% probandCNVs.CG$Inheritance &
                                                                                     "Maternal" %in% probandCNVs.CG$Inheritance)
  
  probandSNVs.CG <- data.frame(CGdata[[i]]$SNVs)
  fatherCNVs.CG <- data.frame(CGdata[[i]]$Father.CNV)
  fatherSNVs.CG <- data.frame(CGdata[[i]]$Father.SNV)
  motherCNVs.CG <- data.frame(CGdata[[i]]$Mother.CNV)
  motherSNVs.CG <- data.frame(CGdata[[i]]$Mother.SNV)
  
  probandSNVs.CG <- probandSNVs.CG[which(probandSNVs.CG$freq_max < 0.1), ]
  fatherSNVs.CG <- fatherSNVs.CG[which(fatherSNVs.CG$freq_max < 0.1), ]
  motherSNVs.CG <- motherSNVs.CG[which(motherSNVs.CG$freq_max < 0.1), ]
  
  if(!("Paternal" %in% probandCNVs.CG$Inheritance & "Maternal" %in% probandCNVs.CG$Inheritance)){
    if(sum(probandCNVs.CG$Inheritance == "Maternal") > 0){
      maternalCount.CG = maternalCount.CG + 1
      #motherCNVs.CG <- GRanges(motherCNVs.CG$chrAnn, IRanges(motherCNVs.CG$STARTAnn, motherCNVs.CG$ENDAnn), "*")
      proband.motherCNVs.CG <- GRanges(probandCNVs.CG$chrAnn[probandCNVs.CG$Inheritance == "Maternal"],
                                     IRanges(probandCNVs.CG$STARTAnn[probandCNVs.CG$Inheritance == "Maternal"],
                                             probandCNVs.CG$ENDAnn[probandCNVs.CG$Inheritance == "Maternal"]), "*")
      motherCNVs.CG <- proband.motherCNVs.CG
    
      if(nrow(probandSNVs.CG) > 0){
        names(probandSNVs.CG)[11] <- "SampleData"
      
        probandSNVs.CG.g <- GRanges(probandSNVs.CG$CHROM, IRanges(probandSNVs.CG$start, probandSNVs.CG$end), "*")
        olap <- findOverlaps(probandSNVs.CG.g, proband.motherCNVs.CG)
        probandSNVs.CG.in.maternalCNVs <- rbind(probandSNVs.CG.in.maternalCNVs, probandSNVs.CG[unique(olap@from), ])
      }
    
      if(nrow(motherSNVs.CG) > 0){
        names(motherSNVs.CG)[11] <- "SampleData"
        motherSNVs.CG <- motherSNVs.CG[,names(motherSNVs.CG)[names(motherSNVs.CG) != "Inheritance.SNV"]]
        motherSNVs.CG.g <- GRanges(motherSNVs.CG$CHROM, IRanges(motherSNVs.CG$start, motherSNVs.CG$end), "*")
        olap <- findOverlaps(motherSNVs.CG.g, motherCNVs.CG)
        motherSNVs.CG.in.maternalCNVs <- rbind(motherSNVs.CG.in.maternalCNVs, motherSNVs.CG[unique(olap@from), ])
      }
    }
  
    if(sum(probandCNVs.CG$Inheritance == "Paternal") > 0){
      paternalCount.CG = paternalCount.CG + 1
      #fatherCNVs.CG <- GRanges(fatherCNVs.CG$chrAnn, IRanges(fatherCNVs.CG$STARTAnn, fatherCNVs.CG$ENDAnn), "*")
      proband.fatherCNVs.CG <- GRanges(probandCNVs.CG$chrAnn[probandCNVs.CG$Inheritance == "Paternal"],
                                     IRanges(probandCNVs.CG$STARTAnn[probandCNVs.CG$Inheritance == "Paternal"],
                                             probandCNVs.CG$ENDAnn[probandCNVs.CG$Inheritance == "Paternal"]), "*")
      fatherCNVs.CG <- proband.fatherCNVs.CG
    
      if(nrow(probandSNVs.CG) > 0){
        names(probandSNVs.CG)[11] <- "SampleData"
      
        probandSNVs.CG.g <- GRanges(probandSNVs.CG$CHROM, IRanges(probandSNVs.CG$start, probandSNVs.CG$end), "*")
        olap <- findOverlaps(probandSNVs.CG.g, proband.fatherCNVs.CG)
        probandSNVs.CG.in.paternalCNVs <- rbind(probandSNVs.CG.in.paternalCNVs, probandSNVs.CG[unique(olap@from), ])
      }
    
      if(nrow(fatherSNVs.CG) > 0){
        names(fatherSNVs.CG)[11] <- "SampleData"
        fatherSNVs.CG <- fatherSNVs.CG[,names(fatherSNVs.CG)[names(fatherSNVs.CG) != "Inheritance.SNV"]]
      
        fatherSNVs.CG.g <- GRanges(fatherSNVs.CG$CHROM, IRanges(fatherSNVs.CG$start, fatherSNVs.CG$end), "*")
        olap <- findOverlaps(fatherSNVs.CG.g, fatherCNVs.CG)
        fatherSNVs.CG.in.paternalCNVs <- rbind(fatherSNVs.CG.in.paternalCNVs, fatherSNVs.CG[unique(olap@from), ])
      }
    }
  }
  message(i)
}

#Write SNVs to tsv
write.table(probandSNVs.CG.in.paternalCNVs, "CG_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(probandSNVs.CG.in.maternalCNVs, "CG_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(fatherSNVs.CG.in.paternalCNVs, "CG_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(motherSNVs.CG.in.maternalCNVs, "CG_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)


##### SSC DATA ##### (1510 CH events)

setwd("C:/Users/acale/Desktop/Fall Co-op 2021/CHASE Files/SSC")

SSCdata <- yaml::yaml.load_file("SSC_CH_Data10P.yaml")

probandSNVs.SSC.in.maternalCNVs <- data.frame()
probandSNVs.SSC.in.paternalCNVs <- data.frame()
motherSNVs.SSC.in.maternalCNVs <- data.frame()
fatherSNVs.SSC.in.paternalCNVs <- data.frame()

maternalCount.SSC = 0
paternalCount.SSC = 0

motherCount.SSC = 0
fatherCount.SSC = 0
probandBothPaternalmaternalCount.SSC = 0

for(i in 1:length(SSCdata)){
  probandCNVs.SSC <- data.frame(SSCdata[[i]]$CNVs)
  
  motherCount.SSC = motherCount.SSC + as.numeric("Maternal" %in% probandCNVs.SSC$Inheritance)
  fatherCount.SSC = fatherCount.SSC + as.numeric("Paternal" %in% probandCNVs.SSC$Inheritance)
  probandBothPaternalmaternalCount.SSC = probandBothPaternalmaternalCount.SSC + as.numeric("Paternal" %in% probandCNVs.SSC$Inheritance &
                                                                                             "Maternal" %in% probandCNVs.SSC$Inheritance)
  
  probandSNVs.SSC <- data.frame(SSCdata[[i]]$SNVs)
  fatherCNVs.SSC <- data.frame(SSCdata[[i]]$Father.CNV)
  fatherSNVs.SSC <- data.frame(SSCdata[[i]]$Father.SNV)
  motherCNVs.SSC <- data.frame(SSCdata[[i]]$Mother.CNV)
  motherSNVs.SSC <- data.frame(SSCdata[[i]]$Mother.SNV)
  
  probandSNVs.SSC <- probandSNVs.SSC[which(probandSNVs.SSC$freq_max < 0.1), ]
  fatherSNVs.SSC <- fatherSNVs.SSC[which(fatherSNVs.SSC$freq_max < 0.1), ]
  motherSNVs.SSC <- motherSNVs.SSC[which(motherSNVs.SSC$freq_max < 0.1), ]
  
  if(!("Paternal" %in% probandCNVs.SSC$Inheritance & "Maternal" %in% probandCNVs.SSC$Inheritance)){
    if(sum(probandCNVs.SSC$Inheritance == "Maternal") > 0){
      maternalCount.SSC = maternalCount.SSC + 1
      #motherCNVs.SSC <- GRanges(motherCNVs.SSC$chrAnn, IRanges(motherCNVs.SSC$STARTAnn, motherCNVs.SSC$ENDAnn), "*")
      proband.motherCNVs.SSC <- GRanges(probandCNVs.SSC$chrAnn[probandCNVs.SSC$Inheritance == "Maternal"],
                                        IRanges(probandCNVs.SSC$STARTAnn[probandCNVs.SSC$Inheritance == "Maternal"],
                                                probandCNVs.SSC$ENDAnn[probandCNVs.SSC$Inheritance == "Maternal"]), "*")
      motherCNVs.SSC <- proband.motherCNVs.SSC
      
      if(nrow(probandSNVs.SSC) > 0){
        names(probandSNVs.SSC)[11] <- "SampleData"
        
        probandSNVs.SSC.g <- GRanges(probandSNVs.SSC$CHROM, IRanges(probandSNVs.SSC$start, probandSNVs.SSC$end), "*")
        olap <- findOverlaps(probandSNVs.SSC.g, proband.motherCNVs.SSC)
        probandSNVs.SSC.in.maternalCNVs <- rbind(probandSNVs.SSC.in.maternalCNVs, probandSNVs.SSC[unique(olap@from), ])
      }
      
      if(nrow(motherSNVs.SSC) > 0){
        names(motherSNVs.SSC)[11] <- "SampleData"
        motherSNVs.SSC <- motherSNVs.SSC[,names(motherSNVs.SSC)[names(motherSNVs.SSC) != "Inheritance.SNV"]]
        motherSNVs.SSC.g <- GRanges(motherSNVs.SSC$CHROM, IRanges(motherSNVs.SSC$start, motherSNVs.SSC$end), "*")
        olap <- findOverlaps(motherSNVs.SSC.g, motherCNVs.SSC)
        motherSNVs.SSC.in.maternalCNVs <- rbind(motherSNVs.SSC.in.maternalCNVs, motherSNVs.SSC[unique(olap@from), ])
      }
    }
    
    if(sum(probandCNVs.SSC$Inheritance == "Paternal") > 0){
      paternalCount.SSC = paternalCount.SSC + 1
      #fatherCNVs.SSC <- GRanges(fatherCNVs.SSC$chrAnn, IRanges(fatherCNVs.SSC$STARTAnn, fatherCNVs.SSC$ENDAnn), "*")
      proband.fatherCNVs.SSC <- GRanges(probandCNVs.SSC$chrAnn[probandCNVs.SSC$Inheritance == "Paternal"],
                                        IRanges(probandCNVs.SSC$STARTAnn[probandCNVs.SSC$Inheritance == "Paternal"],
                                                probandCNVs.SSC$ENDAnn[probandCNVs.SSC$Inheritance == "Paternal"]), "*")
      fatherCNVs.SSC <- proband.fatherCNVs.SSC
      
      if(nrow(probandSNVs.SSC) > 0){
        names(probandSNVs.SSC)[11] <- "SampleData"
        
        probandSNVs.SSC.g <- GRanges(probandSNVs.SSC$CHROM, IRanges(probandSNVs.SSC$start, probandSNVs.SSC$end), "*")
        olap <- findOverlaps(probandSNVs.SSC.g, proband.fatherCNVs.SSC)
        probandSNVs.SSC.in.paternalCNVs <- rbind(probandSNVs.SSC.in.paternalCNVs, probandSNVs.SSC[unique(olap@from), ])
      }
      
      if(nrow(fatherSNVs.SSC) > 0){
        names(fatherSNVs.SSC)[11] <- "SampleData"
        fatherSNVs.SSC <- fatherSNVs.SSC[,names(fatherSNVs.SSC)[names(fatherSNVs.SSC) != "Inheritance.SNV"]]
        
        fatherSNVs.SSC.g <- GRanges(fatherSNVs.SSC$CHROM, IRanges(fatherSNVs.SSC$start, fatherSNVs.SSC$end), "*")
        olap <- findOverlaps(fatherSNVs.SSC.g, fatherCNVs.SSC)
        fatherSNVs.SSC.in.paternalCNVs <- rbind(fatherSNVs.SSC.in.paternalCNVs, fatherSNVs.SSC[unique(olap@from), ])
      }
    }
  }
  message(i)
}

#Write SNVs to tsv
write.table(probandSNVs.SSC.in.paternalCNVs, "SSC_ProbandSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(probandSNVs.SSC.in.maternalCNVs, "SSC_ProbandSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(fatherSNVs.SSC.in.paternalCNVs, "SSC_FatherSNVs_in_PaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(motherSNVs.SSC.in.maternalCNVs, "SSC_MotherSNVs_in_MaternalCNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)

###########################################################################################################################################################

#SNV COUNT BY VARIANT TYPE

##### ILMN DATA ######

#Read SNVs tsv
probandCH_SNVs.paternal <- read.delim("ILMN_ProbandSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F)
probandCH_SNVs.maternal <- read.delim("ILMN_ProbandSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F)
paternalCH_SNVs <- read.delim("ILMN_FatherSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F)
maternalCH_SNVs <- read.delim("ILMN_MotherSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F)

#All variants
probandPaternal_del.ALL <- length(probandCH_SNVs.paternal$Sample.ID[!duplicated(paste(probandCH_SNVs.paternal$X.Sample, probandCH_SNVs.paternal$X.id))])                
probandPaternal_del.ALL.COUNT <- length(unique(probandCH_SNVs.paternal$Sample.ID))   #25 probands
probandMaternal_del.ALL <- length(probandCH_SNVs.maternal$Sample.ID[!duplicated(paste(probandCH_SNVs.maternal$X.Sample, probandCH_SNVs.maternal$X.id))])                
probandMaternal_del.ALL.COUNT <- length(unique(probandCH_SNVs.maternal$Sample.ID))   #26 probands

probandFather.ALL <- length(paternalCH_SNVs$Sample.ID[!duplicated(paste(paternalCH_SNVs$X.Sample, paternalCH_SNVs$X.id))])           
probandFather.ALL.COUNT <- length(unique(paternalCH_SNVs$Sample.ID))                 #15 fathers     
probandMother.ALL <- length(maternalCH_SNVs$Sample.ID[!duplicated(paste(maternalCH_SNVs$X.Sample, maternalCH_SNVs$X.id))])                             
probandMother.ALL.COUNT <- length(unique(maternalCH_SNVs$Sample.ID))                 #20 mothers     


#Synonymous SNVs
synonymousSNVs.paternal <- probandCH_SNVs.paternal[which(probandCH_SNVs.paternal$effect_priority == "synonymous SNV"), ]
synonymousSNVs.maternal <- probandCH_SNVs.maternal[which(probandCH_SNVs.maternal$effect_priority == "synonymous SNV"), ]

probandPaternal_del.SYN <- length(synonymousSNVs.paternal$Sample.ID[!duplicated(paste(synonymousSNVs.paternal$X.Sample, synonymousSNVs.paternal$X.id))])
probandPaternal_del.SYN.COUNT <- length(unique(synonymousSNVs.paternal$Sample.ID))                                                       #10 probands
probandMaternal_del.SYN <- length(synonymousSNVs.maternal$Sample.ID[!duplicated(paste(synonymousSNVs.maternal$X.Sample, synonymousSNVs.maternal$X.id))])
probandMaternal_del.SYN.COUNT <- length(unique(synonymousSNVs.maternal$Sample.ID))                                                       #14 probands

synonymousSNVs.father <- nrow(unique(paternalCH_SNVs[which(paternalCH_SNVs$effect_priority == "synonymous SNV"), 
                                                        c("X.Sample", "X.id")]))           
synonymousSNVs.fatherCOUNT <- length(unique(paternalCH_SNVs$Sample.ID[which(paternalCH_SNVs$effect_priority == "synonymous SNV")]))      #1 fathers
synonymousSNVs.mother <- nrow(unique(maternalCH_SNVs[which(maternalCH_SNVs$effect_priority == "synonymous SNV"), 
                                                     c("X.Sample", "X.id")]))    
synonymousSNVs.motherCOUNT <- length(maternalCH_SNVs$Sample.ID[which(maternalCH_SNVs$effect_priority == "synonymous SNV")])              #5 mothers
 

#Missense SNVs
missenseSNVs.paternal <- probandCH_SNVs.paternal[which(probandCH_SNVs.paternal$effect_priority == "nonsynonymous SNV"), ]
missenseSNVs.maternal <- probandCH_SNVs.maternal[which(probandCH_SNVs.maternal$effect_priority == "nonsynonymous SNV"), ]

probandPaternal_del.MIS <- length(missenseSNVs.paternal$Sample.ID[!duplicated(paste(missenseSNVs.paternal$X.Sample, missenseSNVs.paternal$X.id))])                                                                      
probandPaternal_del.MIS.COUNT <- length(unique(missenseSNVs.paternal$Sample.ID))                                                         #15 probands
probandMaternal_del.MIS <- length(missenseSNVs.maternal$Sample.ID[!duplicated(paste(missenseSNVs.maternal$X.Sample, missenseSNVs.maternal$X.id))])                                                                      
probandMaternal_del.MIS.COUNT <- length(unique(missenseSNVs.maternal$Sample.ID))                                                         #14 probands

missenseSNVs.father <- nrow(unique(paternalCH_SNVs[which(paternalCH_SNVs$effect_priority == "nonsynonymous SNV"), 
                                                   c("X.Sample", "X.id")]))
missenseSNVs.fatherCOUNT <- length(unique(paternalCH_SNVs$Sample.ID[which(paternalCH_SNVs$effect_priority == "nonsynonymous SNV")]))     #11 fathers
missenseSNVs.mother <- nrow(unique(maternalCH_SNVs[which(maternalCH_SNVs$effect_priority == "nonsynonymous SNV"), 
                                                   c("X.Sample", "X.id")]))
missenseSNVs.motherCOUNT <- length(unique(maternalCH_SNVs$Sample.ID[which(maternalCH_SNVs$effect_priority == "nonsynonymous SNV")]))     #17 mothers


#Loss of Function SNVs
LoFSNVs.paternal <- probandCH_SNVs.paternal[which(probandCH_SNVs.paternal$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                           probandCH_SNVs.paternal$typeseq_priority == "splicing"), ]
LoFSNVs.maternal <- probandCH_SNVs.maternal[which(probandCH_SNVs.maternal$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                    probandCH_SNVs.maternal$typeseq_priority == "splicing"), ]

probandPaternal_del.LOF <- length(LoFSNVs.paternal$Sample.ID[!duplicated(paste(LoFSNVs.paternal$X.Sample, LoFSNVs.paternal$X.id))])                                                                            
probandPaternal_del.LOF.COUNT <- length(unique(LoFSNVs.paternal$Sample.ID))                                                              #4 probands
probandMaternal_del.LOF <- length(LoFSNVs.maternal$Sample.ID[!duplicated(paste(LoFSNVs.maternal$X.Sample, LoFSNVs.maternal$X.id))])                                                                          
probandMaternal_del.LOF.COUNT <- length(unique(LoFSNVs.maternal$Sample.ID))                                                              #2 probands

LoFSNVs.father <- nrow(unique(paternalCH_SNVs[which(paternalCH_SNVs$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                      paternalCH_SNVs$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
LoFSNVs.fatherCOUNT <- length(unique(paternalCH_SNVs$Sample.ID[which(paternalCH_SNVs$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                                    paternalCH_SNVs$typeseq_priority == "splicing")]))
LoFSNVs.mother <- nrow(unique(maternalCH_SNVs[which(maternalCH_SNVs$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                      maternalCH_SNVs$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
LoFSNVs.motherCOUNT <- length(unique(maternalCH_SNVs$Sample.ID[!duplicated(paste(maternalCH_SNVs$X.Sample, maternalCH_SNVs$X.id))]
                                     [which(maternalCH_SNVs$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                                    maternalCH_SNVs$typeseq_priority == "splicing")]))

SNVcomparison <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                      proband_paternal_del = c(probandPaternal_del.ALL, probandPaternal_del.SYN, probandPaternal_del.MIS, probandPaternal_del.LOF),
                      proband_maternal_del = c(probandMaternal_del.ALL, probandMaternal_del.SYN, probandMaternal_del.MIS, probandMaternal_del.LOF), 
                      father = c(probandFather.ALL, synonymousSNVs.father, missenseSNVs.father, LoFSNVs.father),
                      mother = c(probandMother.ALL, synonymousSNVs.mother, missenseSNVs.mother, LoFSNVs.mother))

SNVcomparison.df <- as.data.frame(SNVcomparison)

##### CG DATA ######

#Read SNVs tsv
probandCH_SNVs.paternal.CG <- read.delim("CG_ProbandSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F)
probandCH_SNVs.maternal.CG <- read.delim("CG_ProbandSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F)
paternalCH_SNVs.CG <- read.delim("CG_FatherSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F)
maternalCH_SNVs.CG <- read.delim("CG_MotherSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F)

#All variants
probandPaternal_del.ALLCG <- length(probandCH_SNVs.paternal.CG$Sample.ID
                                    [!duplicated(paste(probandCH_SNVs.paternal.CG$X.Sample, probandCH_SNVs.paternal.CG$X.id))])                 
probandPaternal_del.ALL.COUNTCG <- length(unique(probandCH_SNVs.paternal.CG$Sample.ID))   #9 probands
probandMaternal_del.ALLCG <- length(probandCH_SNVs.maternal.CG$Sample.ID
                                    [!duplicated(paste(probandCH_SNVs.maternal.CG$X.Sample, probandCH_SNVs.maternal.CG$X.id))])                 
probandMaternal_del.ALL.COUNTCG <- length(unique(probandCH_SNVs.maternal.CG$Sample.ID))   #5 probands

probandFather.ALLCG <- length(paternalCH_SNVs.CG$Sample.ID
                              [!duplicated(paste(paternalCH_SNVs.CG$X.Sample, paternalCH_SNVs.CG$X.id))])                             
probandFather.ALLCG.COUNT <- length(unique(paternalCH_SNVs.CG$Sample.ID))                 #6 fathers     
probandMother.ALLCG <- length(maternalCH_SNVs.CG$Sample.ID
                              [!duplicated(paste(maternalCH_SNVs.CG$X.Sample, maternalCH_SNVs.CG$X.id))])                               
probandMother.ALLCG.COUNT <- length(unique(maternalCH_SNVs.CG$Sample.ID))                 #2 mothers     


#Synonymous SNVs
synonymousSNVs.paternal.CG <- probandCH_SNVs.paternal.CG[which(probandCH_SNVs.paternal.CG$effect_priority == "synonymous SNV"), ]
synonymousSNVs.maternal.CG <- probandCH_SNVs.maternal.CG[which(probandCH_SNVs.maternal.CG$effect_priority == "synonymous SNV"), ]

probandPaternal_del.SYNCG <- length(synonymousSNVs.paternal.CG$Sample.ID
                                    [!duplicated(paste(synonymousSNVs.paternal.CG$X.Sample, synonymousSNVs.paternal.CG$X.id))])                                                                     
probandPaternal_del.SYNCG.COUNT <- length(unique(synonymousSNVs.paternal.CG$Sample.ID))                                                       #2 probands
probandMaternal_del.SYNCG <- length(synonymousSNVs.maternal.CG$Sample.ID
                                    [!duplicated(paste(synonymousSNVs.maternal.CG$X.Sample, synonymousSNVs.maternal.CG$X.id))])                                                                    
probandMaternal_del.SYNCG.COUNT <- length(unique(synonymousSNVs.maternal.CG$Sample.ID))                                                       #3 probands

synonymousSNVs.fatherCG <- nrow(unique(paternalCH_SNVs.CG[which(paternalCH_SNVs.CG$effect_priority == "synonymous SNV"), 
                                                          c("X.Sample", "X.id")]))                   
synonymousSNVs.fatherCOUNTCG <- length(unique(paternalCH_SNVs.CG$Sample.ID[which(paternalCH_SNVs.CG$effect_priority == "synonymous SNV")]))      #2 fathers
synonymousSNVs.motherCG <- nrow(unique(maternalCH_SNVs.CG[which(maternalCH_SNVs.CG$effect_priority == "synonymous SNV"), 
                                                          c("X.Sample", "X.id")]))                     
synonymousSNVs.motherCOUNTCG <- length(maternalCH_SNVs.CG$Sample.ID[which(maternalCH_SNVs.CG$effect_priority == "synonymous SNV")])              #2 mothers


#Missense SNVs
missenseSNVs.paternal.CG <- probandCH_SNVs.paternal.CG[which(probandCH_SNVs.paternal.CG$effect_priority == "nonsynonymous SNV"), ]
missenseSNVs.maternal.CG <- probandCH_SNVs.maternal.CG[which(probandCH_SNVs.maternal.CG$effect_priority == "nonsynonymous SNV"), ]

probandPaternal_del.MISCG <- length(missenseSNVs.paternal.CG$Sample.ID
                                    [!duplicated(paste(missenseSNVs.paternal.CG$X.Sample, missenseSNVs.paternal.CG$X.id))])                                                                 
probandPaternal_del.MISCG.COUNT <- length(unique(missenseSNVs.paternal.CG$Sample.ID))                                                          #5 probands
probandMaternal_del.MISCG <- length(missenseSNVs.maternal.CG$Sample.ID
                                    [!duplicated(paste(missenseSNVs.maternal.CG$X.Sample, missenseSNVs.maternal.CG$X.id))])                                                                       
probandMaternal_del.MISCG.COUNT <- length(unique(missenseSNVs.maternal.CG$Sample.ID))                                                          #2 probands

missenseSNVs.fatherCG <- nrow(unique(paternalCH_SNVs.CG[which(paternalCH_SNVs.CG$effect_priority == "nonsynonymous SNV"), 
                                                   c("X.Sample", "X.id")]))
missenseSNVs.fatherCOUNTCG <- length(unique(paternalCH_SNVs.CG$Sample.ID[which(paternalCH_SNVs.CG$effect_priority == "nonsynonymous SNV")]))     
missenseSNVs.motherCG <- nrow(unique(maternalCH_SNVs.CG[which(maternalCH_SNVs.CG$effect_priority == "nonsynonymous SNV"), 
                                                        c("X.Sample", "X.id")]))                
missenseSNVs.motherCOUNTCG <- length(unique(maternalCH_SNVs.CG$Sample.ID[which(maternalCH_SNVs.CG$effect_priority == "nonsynonymous SNV")]))     


#Loss of Function SNVs
LoFSNVs.paternal.CG <- probandCH_SNVs.paternal.CG[which(probandCH_SNVs.paternal.CG$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                          probandCH_SNVs.paternal.CG$typeseq_priority == "splicing"), ]
LoFSNVs.maternal.CG <- probandCH_SNVs.maternal.CG[which(probandCH_SNVs.maternal.CG$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                          probandCH_SNVs.maternal.CG$typeseq_priority == "splicing"), ]

probandPaternal_del.LOFCG <- length(LoFSNVs.paternal.CG$Sample.ID
                                    [!duplicated(paste(LoFSNVs.paternal.CG$X.Sample, LoFSNVs.paternal.CG$X.id))])                                                                           
probandPaternal_del.LOFCG.COUNT <- length(unique(LoFSNVs.paternal.CG$Sample.ID))                                                              #2 probands
probandMaternal_del.LOFCG <- length(LoFSNVs.maternal.CG$Sample.ID
                                    [!duplicated(paste(LoFSNVs.maternal.CG$X.Sample, LoFSNVs.maternal.CG$X.id))])                                                                            
probandMaternal_del.LOFCG.COUNT <- length(unique(LoFSNVs.maternal.CG$Sample.ID))                                                              #0 probands

LoFSNVs.fatherCG <- nrow(unique(paternalCH_SNVs.CG[which(paternalCH_SNVs.CG$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                     paternalCH_SNVs.CG$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
LoFSNVs.fatherCOUNTCG <- length(unique(paternalCH_SNVs.CG$Sample.ID[which(paternalCH_SNVs.CG$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                                            paternalCH_SNVs.CG$typeseq_priority == "splicing")]))
LoFSNVs.motherCG <- nrow(unique(maternalCH_SNVs.CG[which(maternalCH_SNVs.CG$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                           maternalCH_SNVs.CG$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
LoFSNVs.motherCOUNTCG <- length(unique(maternalCH_SNVs.CG$Sample.ID[which(maternalCH_SNVs.CG$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                                            maternalCH_SNVs.CG$typeseq_priority == "splicing")]))

SNVcomparison.CG <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                         proband_paternal_del = c(probandPaternal_del.ALLCG, probandPaternal_del.SYNCG, probandPaternal_del.MISCG, probandPaternal_del.LOFCG),
                         proband_maternal_del = c(probandMaternal_del.ALLCG, probandMaternal_del.SYNCG, probandMaternal_del.MISCG, probandMaternal_del.LOFCG), 
                         father = c(probandFather.ALLCG, synonymousSNVs.fatherCG, missenseSNVs.fatherCG, LoFSNVs.fatherCG),
                         mother = c(probandMother.ALLCG, synonymousSNVs.motherCG, missenseSNVs.motherCG, LoFSNVs.motherCG))

SNVcomparison.dfCG <- as.data.frame(SNVcomparison.CG)

##### SSC DATA ######

setwd("C:/Users/acale/Desktop/Fall Co-op 2021/CHASE Files/SSC")

#Read SNVs tsv
probandCH_SNVs.paternal.SSC <- read.delim("SSC_ProbandSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F)
probandCH_SNVs.maternal.SSC <- read.delim("SSC_ProbandSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F)
paternalCH_SNVs.SSC <- read.delim("SSC_FatherSNVs_in_PaternalCNVs.tsv", stringsAsFactors = F)
maternalCH_SNVs.SSC <- read.delim("SSC_MotherSNVs_in_MaternalCNVs.tsv", stringsAsFactors = F)

#All variants
probandPaternal_del.ALLSSC <- length(probandCH_SNVs.paternal.SSC$Sample.ID)                 #49 SNVs
probandPaternal_del.ALL.COUNTSSC <- length(unique(probandCH_SNVs.paternal.SSC$Sample.ID))   #35 probands
probandMaternal_del.ALLSSC <- length(probandCH_SNVs.maternal.SSC$Sample.ID)                 #46 SNVs
probandMaternal_del.ALL.COUNTSSC <- length(unique(probandCH_SNVs.maternal.SSC$Sample.ID))   #34 probands

probandFather.ALLSSC <- length(paternalCH_SNVs.SSC$Sample.ID)                               #56 SNVs
probandFather.ALLSSC.COUNT <- length(unique(paternalCH_SNVs.SSC$Sample.ID))                 #25 fathers     
probandMother.ALLSSC <- length(maternalCH_SNVs.SSC$Sample.ID)                               #52 SNVs
probandMother.ALLSSC.COUNT <- length(unique(maternalCH_SNVs.SSC$Sample.ID))                 #29 mothers     


#Synonymous SNVs
synonymousSNVs.paternal.SSC <- probandCH_SNVs.paternal.SSC[which(probandCH_SNVs.paternal.SSC$effect_priority == "synonymous SNV"), ]
synonymousSNVs.maternal.SSC <- probandCH_SNVs.maternal.SSC[which(probandCH_SNVs.maternal.SSC$effect_priority == "synonymous SNV"), ]

probandPaternal_del.SYNSSC <- length(synonymousSNVs.paternal.SSC$Sample.ID)                                                                   
probandPaternal_del.SYNSSC.COUNT <- length(unique(synonymousSNVs.paternal.SSC$Sample.ID))                                                       #9 probands
probandMaternal_del.SYNSSV <- length(synonymousSNVs.maternal.SSC$Sample.ID)                                                                    
probandMaternal_del.SYNSSV.COUNT <- length(unique(synonymousSNVs.maternal.SSC$Sample.ID))                                                       #13 probands

synonymousSNVs.fatherSSC <- nrow(unique(paternalCH_SNVs.SSC[which(paternalCH_SNVs.SSC$effect_priority == "synonymous SNV"), 
                                                           c("X.Sample", "X.id")]))                   
synonymousSNVs.fatherCOUNTSSC <- length(unique(paternalCH_SNVs.SSC$Sample.ID[which(paternalCH_SNVs.SSC$effect_priority == "synonymous SNV")]))      #15 fathers
synonymousSNVs.motherSSC <- nrow(unique(maternalCH_SNVs.SSC[which(maternalCH_SNVs.SSC$effect_priority == "synonymous SNV"), 
                                                            c("X.Sample", "X.id")]))                  
synonymousSNVs.motherCOUNTSSC <- length(maternalCH_SNVs.SSC$Sample.ID[which(maternalCH_SNVs.SSC$effect_priority == "synonymous SNV")])              #12 mothers


#Missense SNVs
missenseSNVs.paternal.SSC <- probandCH_SNVs.paternal.SSC[which(probandCH_SNVs.paternal.SSC$effect_priority == "nonsynonymous SNV"), ]
missenseSNVs.maternal.SSC <- probandCH_SNVs.maternal.SSC[which(probandCH_SNVs.maternal.SSC$effect_priority == "nonsynonymous SNV"), ]

probandPaternal_del.MISSSC <- length(missenseSNVs.paternal.SSC$Sample.ID)                                                                        
probandPaternal_del.MISSSC.COUNT <- length(unique(missenseSNVs.paternal.SSC$Sample.ID))                                                          #26 probands
probandMaternal_del.MISSSC <- length(missenseSNVs.maternal.SSC$Sample.ID)                                                                        
probandMaternal_del.MISSSC.COUNT <- length(unique(missenseSNVs.maternal.SSC$Sample.ID))                                                          #21 probands

missenseSNVs.fatherSSC <- nrow(unique(paternalCH_SNVs.SSC[which(paternalCH_SNVs.SSC$effect_priority == "nonsynonymous SNV"), 
                                                         c("X.Sample", "X.id")]))                
missenseSNVs.fatherCOUNTSSC <- length(unique(paternalCH_SNVs.SSC$Sample.ID[which(paternalCH_SNVs.SSC$effect_priority == "nonsynonymous SNV")]))     #15 fathers
missenseSNVs.motherSSC <- nrow(unique(maternalCH_SNVs.SSC[which(maternalCH_SNVs.SSC$effect_priority == "nonsynonymous SNV"), 
                                                          c("X.Sample", "X.id")]))                 
missenseSNVs.motherCOUNTSSC <- length(unique(maternalCH_SNVs.SSC$Sample.ID[which(maternalCH_SNVs.SSC$effect_priority == "nonsynonymous SNV")]))     #20 mothers


#Loss of Function SNVs
LoFSNVs.paternal.SSC <- probandCH_SNVs.paternal.SSC[which(probandCH_SNVs.paternal.SSC$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                            probandCH_SNVs.paternal.SSC$typeseq_priority == "splicing"), ]
LoFSNVs.maternal.SSC <- probandCH_SNVs.maternal.SSC[which(probandCH_SNVs.maternal.SSC$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                            probandCH_SNVs.maternal.SSC$typeseq_priority == "splicing"), ]

probandPaternal_del.LOFSSC <- length(LoFSNVs.paternal.SSC$Sample.ID)                                                                            #2 SNVs
probandPaternal_del.LOFSSC.COUNT <- length(unique(LoFSNVs.paternal.SSC$Sample.ID))                                                              #2 probands
probandMaternal_del.LOFSSC <- length(LoFSNVs.maternal.SSC$Sample.ID)                                                                            #4 SNVs
probandMaternal_del.LOFSSC.COUNT <- length(unique(LoFSNVs.maternal.SSC$Sample.ID))                                                              #4 probands

LoFSNVs.fatherSSC <- nrow(unique(paternalCH_SNVs.SSC[which(paternalCH_SNVs.SSC$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                             paternalCH_SNVs.SSC$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
LoFSNVs.fatherCOUNTSSC <- length(unique(paternalCH_SNVs.SSC$Sample.ID[which(paternalCH_SNVs.SSC$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                                              paternalCH_SNVs.SSC$typeseq_priority == "splicing")]))
LoFSNVs.motherSSC <- nrow(unique(maternalCH_SNVs.SSC[which(maternalCH_SNVs.SSC$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                             maternalCH_SNVs.SSC$typeseq_priority == "splicing"), c("X.Sample", "X.id")]))
LoFSNVs.motherCOUNTSSC <- length(unique(maternalCH_SNVs.SSC$Sample.ID[which(maternalCH_SNVs.SSC$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | 
                                                                              maternalCH_SNVs.SSC$typeseq_priority == "splicing")]))

SNVcomparison.SSC <- list(VariantType = c("All Variants", "Synonymous", "Missense", "LoF"),
                          proband_paternal_del = c(probandPaternal_del.ALLSSC, probandPaternal_del.SYNSSC, probandPaternal_del.MISSSC, probandPaternal_del.LOFSSC),
                          proband_maternal_del = c(probandMaternal_del.ALLSSC, probandMaternal_del.SYNSSV, probandMaternal_del.MISSSC, probandMaternal_del.LOFSSC), 
                          father = c(probandFather.ALLSSC, synonymousSNVs.fatherSSC, missenseSNVs.fatherSSC, LoFSNVs.fatherSSC),
                          mother = c(probandMother.ALLSSC, synonymousSNVs.motherSSC, missenseSNVs.motherSSC, LoFSNVs.motherSSC))

SNVcomparison.dfSSC <- as.data.frame(SNVcomparison.SSC)
