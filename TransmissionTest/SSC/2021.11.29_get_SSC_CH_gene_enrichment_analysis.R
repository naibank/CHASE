library(data.table)
meta <- read.delim("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_metadata.tsv", stringsAsFactors = F)
meta$Relation <- gsub("other sibling", "unaffected sibling", meta$Relation)
meta$Relation <- gsub("father|mother", "parent", meta$Relation)

lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

## get list of all_genes
all_genes <- unique(unlist(gsMain))

for (relation in c("proband", "unaffected_sibling")){
  if (relation == "proband"){
    CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_proband_CNV_SNV_table_eventfreq.tsv", data.table = F)
    child = "proband" 
  }
  if (relation == "unaffected_sibling"){
    CNV_SNV_table <- fread("/Users/shaniawu/SickKids CHASE/SSC/data/SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv", data.table = F)
    child = "unaffected sibling"
  }
  
  CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                        CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF"
  
  ###Bank - only use CH events where genes are in gsMain (all_genes)
  table <- CNV_SNV_table[which(CNV_SNV_table$entrez_id %in% all_genes), ]
  
  test.out <- data.frame()
  
  for (gsi in 1:length(gsMain)){
    for(freq in c(0.01, 0.001, 0.0001, 0.00001)){ # event freq cut-offs; 0.01 = no cut-off
      for(pRecx in c(0, 0.5, 0.9)){ # pRec cut-offs
        for(variant_type in c("LoF", "nonsynonymous SNV")){
          if (pRecx == 0){
            #Bank, get other genes for this frequency cut-off and variant type
            other_genes <- all_genes[which(!all_genes %in% table$entrez_id[which(table$event_freq < freq &
                                                                                   table$effect_priority == variant_type)])]
            
            #count of CH genes in gsMain[[gsi]] w/ event_freq and variant_type
            target.CH.ix <- which(table$event_freq < freq &
                                    table$effect_priority == variant_type &
                                    table$entrez_id %in% gsMain[[gsi]])
            target.CH <- length(unique(table$entrez_id[target.CH.ix]))
            
            
            #count of other_genes in gsMain[[gsi]]
            target.other <- length(other_genes[which(other_genes %in% gsMain[[gsi]])])
            
            #count of CH genes not in gsMain[[gsi]] w/ event_freq and variant_type
            bg.CH.ix <- which(table$event_freq < freq &
                                table$effect_priority == variant_type &
                                !table$entrez_id %in% gsMain[[gsi]])
            bg.CH <- length(unique(table$entrez_id[bg.CH.ix]))
            
            #count of other_genes not in gsMain[[gsi]]
            bg.other <- length(other_genes[which(!other_genes %in% gsMain[[gsi]])])
          } else {
            #Bank, get other genes for this frequency cut-off and variant type
            other_genes <- all_genes[which(!all_genes %in% table$entrez_id[which(table$event_freq < freq &
                                                                                   table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                                                                   table$effect_priority == variant_type)])]
            #count of CH genes in gsMain[[gsi]] w/ event_freq and variant_type
            target.CH.ix <- which(table$event_freq < freq &
                                    table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                    table$effect_priority == variant_type &
                                    table$entrez_id %in% gsMain[[gsi]])
            target.CH <- length(unique(table$entrez_id[target.CH.ix]))
            
            #count of other_genes in gsMain[[gsi]]
            target.other <- length(other_genes[which(other_genes %in% gsMain[[gsi]])])
            
            #count of CH genes not in gsMain[[gsi]] w/ event_freq and variant_type
            bg.CH.ix <- which(table$event_freq < freq &
                                table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                                table$effect_priority == variant_type &
                                !table$entrez_id %in% gsMain[[gsi]])
            bg.CH <- length(unique(table$entrez_id[bg.CH.ix]))
            
            #count of other_genes not in gsMain[[gsi]]
            bg.other <- length(other_genes[which(!other_genes %in% gsMain[[gsi]])])
          }
          
          df.test <- data.frame("target" = c(target.CH, target.other),
                                "background" = c(bg.CH, bg.other))
          
          test <- fisher.test(df.test, alternative = "greater")
          test.out <- rbind(test.out, 
                            data.frame("gene_set" = names(gsMain[gsi]), # gene set name
                                       "event_freq" = format(freq, scientific = F),
                                       "pRec" = pRecx,
                                       "variant_type" = variant_type,
                                       "OR" = signif(test$estimate, digits = 3),
                                       "P" = signif(test$p.value, digits = 3),
                                       target.CH,
                                       target.other,
                                       bg.CH,
                                       bg.other))
        }
      }
    }
  }
  table_name <- sprintf("SSC_CH_%s_gene_enrichment_analysis_gsMain.tsv", relation)
  write.table(test.out, table_name, sep="\t", row.names=F, quote=F, col.names=T)
}
