library(data.table)
library(survival)

source("functions.R")

probandTable <- fread("SNV_in_CH_tables/Proband_parents.tsv", data.table = F)
unaffTable <- fread("SNV_in_CH_tables/Unaff_parents.tsv", data.table = F)

probandTable <- probandTable[which(probandTable$gnomAD_oe_lof_upper > 0.35 | is.na(probandTable$gnomAD_oe_lof_upper)), ]
unaffTable <- unaffTable[which(unaffTable$gnomAD_oe_lof_upper > 0.35 | is.na(unaffTable$gnomAD_oe_lof_upper)), ]

ge <- readLines("../recessive_genes/GE_Neurodevelopment_biallelic_Xavier_genes.txt")

metatable <- rbind(read.delim("../prerun_family_QCs/MSSNG_CG.metadata.tsv"),
                   read.delim("../prerun_family_QCs/MSSNG_ILMN.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_1.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_2.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_3.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_4.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SSC.metadata.tsv"))
trios <- metatable[metatable$family_member == "trios" & metatable$Father.ID != "-" & metatable$Mother.ID != "-", ]


metatable <- rbind(trios, metatable[metatable$Sample.ID %in% c(trios$Father.ID, trios$Mother.ID), ])
affected_parent_fam <- metatable$Family.ID[metatable$Relation %in% c("father", "mother") & metatable$Affection %in% c(0, 2)]

metatable <- metatable[!metatable$Family.ID %in% 
                         affected_parent_fam, ]
probandTable <- probandTable[probandTable$X.Sample %in% metatable$Sample.ID, ]
unaffTable <- unaffTable[unaffTable$X.Sample %in% metatable$Sample.ID, ]

probandTable$effect_priority[which(probandTable$effect_priority == "nonsynonymous SNV")] <- "nonsynonymous SNV"
unaffTable$effect_priority[which(unaffTable$effect_priority == "nonsynonymous SNV")] <- "nonsynonymous SNV"

probandTable$effect_priority[which(probandTable$effect_priority == "nonsynonymous SNV" &
                                     probandTable$damaging_missense_count >= 4)] <- "damaging missense"
unaffTable$effect_priority[which(unaffTable$effect_priority == "nonsynonymous SNV" &
                                   unaffTable$damaging_missense_count >= 4)] <- "damaging missense"

probandTable$effect_priority[which(probandTable$LoF)] <- "LoF"
unaffTable$effect_priority[which(unaffTable$LoF)] <- "LoF"

dt.out <- data.frame()
for(kid in c("proband", "unaffected")){
  if(kid == "proband"){
    subtable <- probandTable
    submeta <- metatable[!metatable$Relation %in% c("father", "mother") & metatable$Affection == 2, ]
    submeta <- rbind(submeta, metatable[metatable$Sample.ID %in% c(submeta$Father.ID, submeta$Mother.ID), ])
  }else{
    subtable <- unaffTable
    submeta <- metatable[!metatable$Relation %in% c("father", "mother") & metatable$Affection == 1, ]
    submeta <- rbind(submeta, metatable[metatable$Sample.ID %in% c(submeta$Father.ID, submeta$Mother.ID), ])
    submeta$Affection[submeta$Relation %in% c("other", "other sibling", "child", "sibling", "unaffected sibling") &
                        submeta$Affection == 1] <- 2
  }
  
for(var in c("synonymous SNV", "nonsynonymous SNV", "damaging missense", "LoF", "damaging variant")){
  # for(freq in c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001)){
  if(var == "nonsynonymous SNV"){
    vartype = c(var, "damaging missense")
  }else if(var == "damaging variant"){
    vartype = c("LoF", "damaging missense")
  }else{
    vartype = var
  }
  
  for(freq in c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005)){
    for(pRec in list(c(0,0.5), c(0.5,0.9), c(0.9,1))){
      
      
      res <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & 
                                         subtable$del_freq_max/100 * subtable$freq_max <= freq & 
                                         subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta)
           
      res.mssng <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & 
                                               subtable$del_freq_max/100 * subtable$freq_max <= freq & 
                                               subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta[submeta$Dataset == "MSSNG", ])
      
      res.ssc <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & 
                                             subtable$del_freq_max/100 * subtable$freq_max <= freq & 
                                               subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta[submeta$Dataset == "SSC", ])
      
      res.spark <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & 
                                               subtable$del_freq_max/100 * subtable$freq_max <= freq & 
                                               subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta[grep("SPARK", submeta$Dataset), ])
      
      names(res.mssng) <- paste("MSSNG", names(res.mssng), sep=".")
      names(res.ssc) <- paste("SSC", names(res.ssc), sep=".")
      names(res.spark) <- paste("SPARK", names(res.spark), sep=".")
      
      if(nrow(res) > 0){
        dt.out <- rbind(dt.out, 
                        data.frame(kid, var, freq, "pRec"=pRec[1], "gene_panel" = F, res, res.mssng, res.ssc, res.spark))
      }             
  
      ###GE gene panel
      if(pRec[1] == 0){
        res <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & subtable$gene_symbol %in% ge &
                                           subtable$del_freq_max/100 * subtable$freq_max <= freq), ], submeta)
        
        res.mssng <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & subtable$gene_symbol %in% ge & 
                                                 subtable$del_freq_max/100 * subtable$freq_max <= freq), ], submeta[submeta$Dataset == "MSSNG", ])
        
        res.ssc <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & subtable$gene_symbol %in% ge & 
                                               subtable$del_freq_max/100 * subtable$freq_max <= freq), ], submeta[submeta$Dataset == "SSC", ])
        
        res.spark <- fisherTest(subtable[which(subtable$effect_priority %in% vartype & subtable$gene_symbol %in% ge & 
                                                 subtable$del_freq_max/100 * subtable$freq_max <= freq), ], submeta[grep("SPARK", submeta$Dataset), ])
        
        names(res.mssng) <- paste("MSSNG", names(res.mssng), sep=".")
        names(res.ssc) <- paste("SSC", names(res.ssc), sep=".")
        names(res.spark) <- paste("SPARK", names(res.spark), sep=".")
        
        if(nrow(res) > 0){
          dt.out <- rbind(dt.out, 
                          data.frame(kid, var, freq, "pRec"=pRec[1], "gene_panel" = T, res, res.mssng, res.ssc, res.spark))
        }      
      }
      ################
    }
  }
}
}

dt.out[,1:8]
write.table(dt.out, "fisher.results.tsv", sep="\t", row.names=F, quote=F, col.names=T)
