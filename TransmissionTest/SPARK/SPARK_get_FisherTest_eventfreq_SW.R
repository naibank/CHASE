library(data.table)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARKAnalysis/")

meta <- rbind(read.delim("./Data/SPARK_WGS_1_metadata_relfixed.tsv", stringsAsFactors = F),
              read.delim("./Data/SPARK_WGS_2_metadata_relfixed.tsv", stringsAsFactors = F),
              read.delim("./Data/SPARK_WGS_3_metadata_relfixed.tsv", stringsAsFactors = F))

for (relation in c("proband", "unaffected_sibling")){
  if (relation == "proband"){
    CNV_SNV_table_full <- fread("./TT/SPARK_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
    
    child = "proband or affected sibling" 
  }
  if (relation == "unaffected_sibling"){
    CNV_SNV_table_full <- fread("./TT/SPARK_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
    
    child = "unaffected sibling"
  }
  for (snv_freq in c(1, 0.1, 0.01)){
    #Adjust SNV freq
    CNV_SNV_table <- CNV_SNV_table_full[CNV_SNV_table_full$freq_max <= snv_freq,]
    #CNV_SNV_table <- merge(CNV_SNV_table, meta[, c("Sample.ID", "Relation")], by.x = "#Sample", by.y = "Sample.ID", all.x = T)
    ## group LOF variants in effect_priority
    CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                          CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank
    
    lof <- read.delim("gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)
    
    test.out <- data.frame()
    
    for(freq in c(0.01, 0.001, 0.005, 0.0005, 0.0001)){ # event freq cut-offs
      for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
        for(variant_type in c("damaging_missense", "LoF", "nonsynonymous","synonymous")){
          
          if (!pRecx == 0){ 
            if (variant_type == "damaging_missense"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                   CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$damaging_missense_count > 3)
            }
            if (variant_type == "LoF"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                   CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="LoF")
            }
            if (variant_type == "nonsynonymous"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                   CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="nonsynonymous SNV")
            }
            if (variant_type == "synonymous"){
              target.ix <- which(CNV_SNV_table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                   CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="synonymous SNV")
            }
            
          }else{ # no pRec cut-off for pRecx == 0
            if (variant_type == "damaging_missense"){
              target.ix <- which(CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$damaging_missense_count > 3)
            }
            if (variant_type == "LoF"){
              target.ix <- which(CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="LoF")
            }
            if (variant_type == "nonsynonymous"){
              target.ix <- which(CNV_SNV_table$event_freq < freq &
                                   CNV_SNV_table$effect_priority =="nonsynonymous SNV")
            }
            if (variant_type == "synonymous"){
              target.ix <- which(CNV_SNV_table$event_freq < freq &
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
                            data.frame("event_freq" = format(freq, scientific = F),
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
    }
    table_name <- paste("SPARK_", relation, sprintf("_FisherTest_eventfreq_snvfreq%s.tsv", snv_freq), sep = "")
    write.table(test.out, paste('./TT/',table_name), sep="\t", row.names=F, quote=F, col.names=T)
    }
}

# #### Find Genes ####
# Find_Gene_Count <- function(full_table, relation, gene) {
#   gene_count <- length(full_table$CHROM[full_table$Relation %in% relation 
#                                         & full_table$gene_symbol %in% gene
#                                         & !full_table$`#id` %in% full_table$`#id`[full_table$Relation == 'proband or affected sibling']])
#   
#   return (gene_count)
# }

####Genes####
# CNV_SNV_table_SSC_proband <- fread("SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
# CNV_SNV_table_SSC_proband <- CNV_SNV_table_SSC_proband[CNV_SNV_table_SSC_proband$freq_max <= 0.01,]
# CNV_SNV_table_SSC_proband$effect_priority[which(CNV_SNV_table_SSC_proband$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
#                                                   CNV_SNV_table_SSC_proband$typeseq_priority == "splicing")] <- "LoF" ###Bank
# 
# CNV_SNV_table_SSC_US <- fread("SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
# CNV_SNV_table_SSC_US <- CNV_SNV_table_SSC_US[CNV_SNV_table_SSC_US$freq_max <= 0.01,]
# CNV_SNV_table_SSC_US$effect_priority[which(CNV_SNV_table_SSC_US$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
#                                              CNV_SNV_table_SSC_US$typeseq_priority == "splicing")] <- "LoF" ###Bank
# 
# 
# SSC_proband_nonsyn_genes_TT <- CNV_SNV_table_SSC_proband[CNV_SNV_table_SSC_proband$effect_priority == 'nonsynonymous SNV' & 
#                                                            CNV_SNV_table_SSC_proband$Relation == 'proband or affected sibling',
#                                                          c('#Sample','gene_symbol')]
# SSC_proband_nonsyn_genes_TT <- data.frame(SSC_proband_nonsyn_genes_TT %>% 
#                                             distinct(.keep_all=T) %>%
#                                             group_by(gene_symbol) %>% 
#                                             filter(n() > 1) %>%
#                                             summarize(Freq=n(), 
#                                                       Parent=Find_Gene_Count(CNV_SNV_table_SSC_proband,'parent',gene_symbol),
#                                                       `Unaffected Sibling`=Find_Gene_Count(CNV_SNV_table_SSC_US,'unaffected sibling',gene_symbol)) %>% 
#                                             arrange(-Freq))
# write.table(SSC_proband_nonsyn_genes_TT, "../TT/Combined Data/SSC/SNV0.01/SSC_proband_nonsyn_genes_TT_snvfreq0.01.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# 
# 
# SSC_proband_lof_genes_TT <- CNV_SNV_table_SSC_proband[CNV_SNV_table_SSC_proband$effect_priority == 'LoF' & 
#                                                         CNV_SNV_table_SSC_proband$Relation == 'proband or affected sibling',
#                                                       c('#Sample','gene_symbol')]
# SSC_proband_lof_genes_TT <- data.frame(SSC_proband_lof_genes_TT %>% 
#                                          distinct(.keep_all=T) %>%
#                                          group_by(gene_symbol) %>% 
#                                          filter(n() > 1) %>%
#                                          summarize(Freq=n(), 
#                                                    Parent=Find_Gene_Count(CNV_SNV_table_SSC_proband,'parent',gene_symbol),
#                                                    `Unaffected Sibling`=Find_Gene_Count(CNV_SNV_table_SSC_US,'unaffected sibling',gene_symbol)) %>% 
#                                          arrange(-Freq))
# write.table(SSC_proband_lof_genes_TT, "../TT/Combined Data/SSC/SNV0.01/SSC_proband_lof_genes_TT_snvfreq0.01.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# 