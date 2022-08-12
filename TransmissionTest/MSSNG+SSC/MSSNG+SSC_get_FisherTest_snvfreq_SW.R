################################################################################################################################################################################
# 
# MSSNG+SSC_get_FisherTest_snvfreq_SW.R
# purpose: outputs Fisher's exact test results for MSSNG+SSC combined parent-proband and 
#           SSC parent-unaffected sibling for different pRec & SNV frequency cut-offs 
#           (includes OR and transmissions rate)
# input: *dataset*_CNV_SNV_table_eventfreq.tsv 
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data
#             -> output of MSSNG+SSC_process_CNV_SNV_table_SW.R
# output: Fisher's exact test results tables (e.g. MSSNG.SSC_parent_proband_FisherTest_snvfreq0.1.tsv)
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data
#             -> used as input for create_TT_results_plot_SW.R
#         
# notes:
#
##############################################################################################################################################################################

library(data.table)
library(dplyr)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data")

lof <- read.delim("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/gene_data/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)


Get_Fisher_Table <- function(CNV_SNV_table_full, name, relation){
  if (relation == "proband"){
    child = "proband or affected sibling" 
  }
  if (relation == "unaffected_sibling"){
    child = "unaffected sibling"
  }
  for (snv_freq in c(1, 0.1, 0.01)){
    #Adjust SNV freq
    CNV_SNV_table <- CNV_SNV_table_full[CNV_SNV_table_full$freq_max <= snv_freq,]

    ## group LOF variants in effect_priority
    CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                          CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank
    
    test.out <- data.frame()
    
    # for(freq in c(0.01, 0.001, 0.005, 0.0005, 0.0001)){ # event freq cut-offs
      for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
        for(variant_type in c("damaging_missense", "LoF", "nonsynonymous","synonymous")){
          
          if (!pRecx == 0){ 
            if (variant_type == "damaging_missense"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                   # CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$damaging_missense_count > 3)
            }
            if (variant_type == "LoF"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                   # CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="LoF")
            }
            if (variant_type == "nonsynonymous"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                   # CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="nonsynonymous SNV")
            }
            if (variant_type == "synonymous"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[which(lof$pRec > pRecx)] &
                                   # CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="synonymous SNV")
            }
            
          }else{ # no pRec cut-off for pRecx == 0
            if (variant_type == "damaging_missense"){
              target.ix <- which(#CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$damaging_missense_count > 3)
            }
            if (variant_type == "LoF"){
              target.ix <- which(#CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="LoF")
            }
            if (variant_type == "nonsynonymous"){
              target.ix <- which(#CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="nonsynonymous SNV")
            }
            if (variant_type == "synonymous"){
              target.ix <- which(#CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="synonymous SNV")
            }
          }
          
          target.child <- nrow(unique(CNV_SNV_table[intersect(target.ix, which(CNV_SNV_table$Relation == child)), 
                                                    c("#Sample", "#id")]))
          #sum(CNV_SNV_table$Relation[target.ix] == child)
          target.parent <- sum(CNV_SNV_table$Relation[target.ix] == "parent")
          
          bg.child <- sum(CNV_SNV_table$Relation[-target.ix] == child)
          bg.parent <- sum(CNV_SNV_table$Relation[-target.ix] == "parent")
          
          df.test <- data.frame("target" = c(target.child, target.parent-target.child),
                                "background" = c(bg.child, bg.parent-bg.child))
          
          test <- fisher.test(df.test, alternative = "greater")
          binom <- binom.test(c(df.test$target[1], df.test$target[2]))
          
          test.two.side <- fisher.test(df.test, alternative = "two.sided")
          test.out <- rbind(test.out, 
                            data.frame(#"event_freq" = format(freq, scientific = F),
                                       "pRec" = pRecx,
                                       "variant_type" = variant_type,
                                       "OR" = signif(test$estimate, digits = 3),
                                       "OR_upper" = signif(test.two.side$conf.int[2], digits = 3),
                                       "OR_lower" = signif(test.two.side$conf.int[1], digits = 3),
                                       "P" = signif(test$p.value, digits = 3),
                                       "rate" = signif(binom$estimate, digits = 2),
                                       "rate_upper" = signif(binom$conf.int[2], digits = 2),
                                       "rate_lower" = signif(binom$conf.int[1], digits = 2),
                                       target.child,
                                       target.parent,
                                       bg.child,
                                       bg.parent))
        }
      }
    #}
    table_name <- sprintf("%s_FisherTest_snvfreq%s.tsv", name, snv_freq)
    write.table(test.out, table_name, sep="\t", row.names=F, quote=F, col.names=T)
  }
}


##############################################################################################################################################################################
### MSSNG+SSC Parent-Proband ####
# meta <- read.delim("../../../data/MSSNG_metadata.tsv", stringsAsFactors = F)
# ssc <- read.delim("../../../data/SSC/SSC_metadata.tsv", stringsAsFactors = F)
# col <- intersect(names(meta), names(ssc))
# 
# meta <- rbind(meta[, col], ssc[, col])
# 
# meta <- meta[, c("Family.ID", "Sample.ID", "Relation")]

MSSNG.SSC_CNV_SNV_table <- data.table::fread("MSSNG.SSC_CNV_SNV_table_eventfreq.tsv", data.table = F)
MSSNG.SSC_CNV_SNV_table <- MSSNG.SSC_CNV_SNV_table[which(MSSNG.SSC_CNV_SNV_table$CNV_freq < 0.01), ]

# bgrates <- data.frame()
# for(snvfreq in c(1, 0.1, 0.01)){
#   data <- MSSNG.SSC_CNV_SNV_table[-grep("SS", MSSNG.SSC_CNV_SNV_table$`#Sample`), ]
#   parent.data <- sum(data$Relation == "parent" & data$freq_max < snvfreq)
#   kid.data <- sum(data$Relation != "parent" & data$freq_max < snvfreq)
#  
#   bgrates <- rbind(bgrates, data.frame(snvfreq, "set" = "proband", "transmission_rate" =  kid.data/parent.data))
# }
# MSSNG.SSC_CNV_SNV_table <- merge(MSSNG.SSC_CNV_SNV_table, meta, by.x = "#Sample", by.y = "Sample.ID", all.x = T)
# 
# dt.out <- data.frame()
# for(kidid in unique(MSSNG.SSC_CNV_SNV_table$sample)){
#   parent <- sum(MSSNG.SSC_CNV_SNV_table$Relation %in% c("parent") & MSSNG.SSC_CNV_SNV_table$sample == kidid)
#   kid <- sum(!MSSNG.SSC_CNV_SNV_table$Relation %in% c("parent") & MSSNG.SSC_CNV_SNV_table$sample == kidid)
# 
#   dt.out <- rbind(dt.out, data.frame(fam, kid/parent))
# }

Get_Fisher_Table(MSSNG.SSC_CNV_SNV_table, name = "MSSNG.SSC_parent_proband", relation = "proband")


##############################################################################################################################################################################
### SSC Parent-Unaffected sibling ####

SSC_CNV_SNV_table <- fread("SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
SSC_CNV_SNV_table <- SSC_CNV_SNV_table[which(SSC_CNV_SNV_table$CNV_freq < 0.01), ]

# for(snvfreq in c(1, 0.1, 0.01)){
#   unaff <- SSC_CNV_SNV_table
#   parent.unaff <- sum(unaff$Relation == "parent" & unaff$freq_max < snvfreq)
#   kid.unaff <- sum(unaff$Relation != "parent" & unaff$freq_max < snvfreq)
#   bgrates <- rbind(bgrates, data.frame(snvfreq, "set" = "SSC_unaffected", "transmission_rate" =  kid.unaff/parent.unaff))
# }

Get_Fisher_Table(SSC_CNV_SNV_table, name = "SSC_parent_US", relation = "unaffected_sibling")

# dt.out.unaff <- data.frame()
# for(kidid in unique(SSC_CNV_SNV_table$sample)){
#   parent <- sum(SSC_CNV_SNV_table$Relation %in% c("parent") & SSC_CNV_SNV_table$sample == kidid)
#   kid <- sum(!SSC_CNV_SNV_table$Relation %in% c("parent") & SSC_CNV_SNV_table$sample == kidid)
#   
#   dt.out.unaff <- rbind(dt.out.unaff, data.frame(fam, kid/parent))
# }

