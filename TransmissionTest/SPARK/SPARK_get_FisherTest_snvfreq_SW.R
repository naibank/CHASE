################################################################################################################################################################################
# 
# SPARK_get_FisherTest_snvfreq_SW.R
# purpose: outputs Fisher's exact test results for SPARK combined parent-proband 
#           for different pRec & SNV frequency cut-offs  
# input: SPARK_CNV_SNV_table_eventfreq.tsv, 
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/SPARK/SPARK_TDT/data
#             -> output of SPARK_process_CNV_SNV_table_SW.R
#        SPARK_WGS_*1-3*_metadata_relfixed.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/Data/
#         
# output: Fisher's exact test results tables (e.g. SPARK_proband_FisherTest_snvfreq1.tsv)
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/SPARK/SPARK_TDT/data
#             -> used as input for SPARK_create_TT_results_plot_SW.R
#         
# notes:
#
##############################################################################################################################################################################

library(data.table)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/")

meta <- rbind(read.delim("./Data/SPARK_WGS_1_metadata_relfixed.tsv", stringsAsFactors = F),
              read.delim("./Data/SPARK_WGS_2_metadata_relfixed.tsv", stringsAsFactors = F),
              read.delim("./Data/SPARK_WGS_3_metadata_relfixed.tsv", stringsAsFactors = F))

for (relation in c("proband", "unaffected_sibling")){
  if (relation == "proband"){
    CNV_SNV_table <- fread("./TT/SPARK_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
    
    child = "proband or affected sibling"
  }
  if (relation == "unaffected_sibling"){
    CNV_SNV_table <- fread("./TT/SPARK_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
    
    child = "unaffected sibling"
  }

  #CNV_SNV_table <- merge(CNV_SNV_table, meta[, c("Sample.ID", "Relation")], by.x = "#Sample", by.y = "Sample.ID", all.x = T)
  ## group LOF variants in effect_priority
  CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                        CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank
  
  lof <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)
  
  for(freq in c(1, 0.1, 0.01)){ # SNV freq ($freq_max) cut-offs 
    test.out <- data.frame()
    #CNV_SNV_table <- CNV_SNV_table_full[CNV_SNV_table_full$freq_max <= snv_freq,]
    for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
      for(variant_type in c("damaging_missense", "LoF", "nonsynonymous","synonymous")){
        
        if (!pRecx == 0){ 
          if (variant_type == "damaging_missense"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                 CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$damaging_missense_count > 3)
          }
          if (variant_type == "LoF"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                 CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$effect_priority =="LoF")
          }
          if (variant_type == "nonsynonymous"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                 CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$effect_priority =="nonsynonymous SNV")
          }
          if (variant_type == "synonymous"){
            target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                 CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$effect_priority =="synonymous SNV")
          }
          
        }else{ # no pRec cut-off for pRecx == 0
          if (variant_type == "damaging_missense"){
            target.ix <- which(CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$damaging_missense_count > 3)
          }
          if (variant_type == "LoF"){
            target.ix <- which(CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$effect_priority =="LoF")
          }
          if (variant_type == "nonsynonymous"){
            target.ix <- which(CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$effect_priority =="nonsynonymous SNV")
          }
          if (variant_type == "synonymous"){
            target.ix <- which(CNV_SNV_table$freq_max <= freq &
                                 CNV_SNV_table$effect_priority =="synonymous SNV")
          }
        }
        
        target.child <- length(unique(CNV_SNV_table[intersect(target.ix, which(CNV_SNV_table$Relation == child)), c("UID")]))
        target.parent <- length(unique(CNV_SNV_table[intersect(target.ix, which(CNV_SNV_table$Relation == "parent")), c("UID")]))
        
        bg.child <- length(unique(CNV_SNV_table[setdiff(which(CNV_SNV_table$Relation == child), target.ix), c("UID")]))
        bg.parent <- length(unique(CNV_SNV_table[setdiff(which(CNV_SNV_table$Relation == "parent"), target.ix), c("UID")]))
        
        df.test <- data.frame("target" = c(target.child, target.parent-target.child),
                              "background" = c(bg.child, bg.parent-bg.child))
        
        test <- fisher.test(df.test, alternative = "greater")
        test.out <- rbind(test.out, 
                          data.frame("snv_freq" = format(freq, scientific = F),
                                     "pRec" = pRecx,
                                     "variant_type" = variant_type,
                                     "OR" = signif(test$estimate, digits = 3),
                                     "P" = signif(test$p.value, digits = 3),
                                     target.child,
                                     target.parent,
                                     bg.child,
                                     bg.parent))
      }
    }
    table_name <- paste("SPARK_", relation, sprintf("_FisherTest_snvfreq%s.tsv", freq), sep = "")
    write.table(test.out, paste('./TT/',table_name), sep="\t", row.names=F, quote=F, col.names=T)
  }
}
