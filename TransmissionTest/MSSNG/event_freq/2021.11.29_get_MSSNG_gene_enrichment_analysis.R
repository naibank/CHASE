library(data.table)
CNV_SNV_table <- data.table::fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/CNV_SNV_table_eventfreq.tsv", data.table = F)

## group LOF variants in effect_priority
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                      CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank

## add "family_type" col indicating whether sample comes from MPX or SPX family
metadata <- fread("/Volumes/T5 ExFAT/F21 Co-op SickKids/CHASE (Bioinformtics Project)/Introduction data/data/MSSNG_metadata.tsv", data.table = F)

mpx <- metadata$`Sample ID`[which(metadata$`Family type` == "MPX")] # MPX sampleIDs 
spx <- metadata$`Sample ID`[which(metadata$`Family type` == "SPX")] # SPX sampleIDs 

CNV_SNV_table$'family_type' <- NA
CNV_SNV_table$family_type[which(CNV_SNV_table$`#Sample` %in% mpx)] <- "MPX"
CNV_SNV_table$family_type[which(CNV_SNV_table$`#Sample` %in% spx)] <- "SPX"
CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$Relation == "proband or affected sibling"), ] ###Bank - only use CH events in proband or affected sibling
#write.table(CNV_SNV_table, "MSSNG_CNV_SNV_table_fam_type.tsv", sep="\t", row.names=F, quote=F, col.names=T)


## get list of all_genes
# CH_genes <- unique(CNV_SNV_table$entrez_id) #some are not in gsMain # Bank - good catch! therefore, we will only use the ones in gsMain
all_genes <- unique(unlist(gsMain))
# other_genes <- all_genes[which(!all_genes %in% CH_genes)] # Bank - other genes will be changed depends on a set of CH events tested, i.e. freq < 0.001


lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)


for (group in c("all", "MPX", "SPX")){ # filter group type
  ###Bank - only use CH events where genes are in gsMain (all_genes)
  if (group == "MPX"){
    table <- CNV_SNV_table[which(CNV_SNV_table$entrez_id %in% all_genes &
                                           CNV_SNV_table$family_type == group), ]
  }
  if (group == "SPX"){
    table <- CNV_SNV_table[which(CNV_SNV_table$entrez_id %in% all_genes &
                                           CNV_SNV_table$family_type == group), ]
  }
  if (group == "all"){
    table <- CNV_SNV_table[which(CNV_SNV_table$entrez_id %in% all_genes), ] 
  }
  
  ## make table
  test.out <- data.frame()
  
  for (gsi in 1:length(gsMain)){
    for(freq in c(0.005, 0.0001, 0.00001)){ # event freq cut-offs; 0.005 = no cut-off
      for(pRecx in c(0, 0.5, 0.9)){ #pRec cut-offs
        for(variant_type in c("LoF", "nonsynonymous SNV")){
          if (!pRecx == 0){ # apply pRec cut-offs
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
            
          } else { # pRecx  == 0, no cut-off
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
  table_name <- sprintf("MSSNG_CH_%s_gene_enrichment_analysis_gsMain.tsv", group)
  write.table(test.out, table_name, sep="\t", row.names=F, quote=F, col.names=T)
}
  

