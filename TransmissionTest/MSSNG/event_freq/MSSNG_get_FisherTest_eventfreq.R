library(data.table)
library(dplyr)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

CNV_SNV_table <- data.table::fread("MSSNG_CNV_SNV_table_eventfreq.tsv", data.table = F)

#Adjust SNV freq
CNV_SNV_table <- CNV_SNV_table[CNV_SNV_table$freq_max <= 0.01,]
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
      
      target.child <- nrow(unique(CNV_SNV_table[intersect(target.ix, which(CNV_SNV_table$Relation == "proband or affected sibling")), 
                                                c("#Sample", "#id")]))
        #sum(CNV_SNV_table$Relation[target.ix] == "proband or affected sibling")
      target.parent <- sum(CNV_SNV_table$Relation[target.ix] == "parent")
      
      bg.child <- sum(CNV_SNV_table$Relation[-target.ix] == "proband or affected sibling")
      bg.parent <- sum(CNV_SNV_table$Relation[-target.ix] == "parent")
      
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

write.table(test.out, "../TT/Combined Data/MSSNG/SNV0.01/MSSNG_FisherTest_eventfreq_snvfreq0.01.tsv", sep="\t", row.names=F, quote=F, col.names=T)

#### Find Genes ####
Find_Gene_Count <- function(full_table, relation, gene) {
  gene_count <- length(full_table$CHROM[full_table$Relation %in% relation 
                                        & full_table$gene_symbol %in% gene
                                        & !full_table$`#id` %in% full_table$`#id`[full_table$Relation == 'proband or affected sibling']])
  
  return (gene_count)
}


MSSNG_proband_nonsyn_genes_TT <- CNV_SNV_table[CNV_SNV_table$effect_priority == 'nonsynonymous SNV' & 
                                                      CNV_SNV_table$Relation == 'proband or affected sibling',
                                               c('#Sample','gene_symbol')]
MSSNG_proband_nonsyn_genes_TT <- data.frame(MSSNG_proband_nonsyn_genes_TT %>% 
                                              distinct(.keep_all=T) %>%
                                              group_by(gene_symbol) %>% 
                                              filter(n() > 1) %>%
                                              summarize(Freq=n(), 
                                                        Parent=Find_Gene_Count(CNV_SNV_table,'parent',gene_symbol)) %>% 
                                              arrange(-Freq))
write.table(MSSNG_proband_nonsyn_genes_TT, "../TT/Combined Data/MSSNG/SNV0.01/MSSNG_proband_nonsyn_genes_TT_snvfreq0.01.tsv", sep="\t", row.names=F, quote=F, col.names=T)

MSSNG_proband_lof_genes_TT <- CNV_SNV_table[CNV_SNV_table$effect_priority == 'LoF' & 
                                                 CNV_SNV_table$Relation == 'proband or affected sibling',
                                               c('#Sample','gene_symbol')]
MSSNG_proband_lof_genes_TT <- data.frame(MSSNG_proband_lof_genes_TT %>% 
                                              distinct(.keep_all=T) %>%
                                              group_by(gene_symbol) %>% 
                                              filter(n() > 1) %>%
                                              summarize(Freq=n(), 
                                                        Parent=Find_Gene_Count(CNV_SNV_table,'parent',gene_symbol)) %>% 
                                              arrange(-Freq))
write.table(MSSNG_proband_lof_genes_TT, "../TT/Combined Data/MSSNG/SNV0.01/MSSNG_proband_lof_genes_TT_snvfreq0.01.tsv", sep="\t", row.names=F, quote=F, col.names=T)


  