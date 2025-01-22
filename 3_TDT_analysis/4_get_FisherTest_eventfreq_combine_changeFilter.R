library(data.table)
library(dplyr)

setwd("CHASE/3_TDT_analysis")

Get_Fisher_Table <- function(CNV_SNV_table_proband_ori, CNV_SNV_table_unaff_sib_ori, GEgenes, name){
  child <- "proband or affected sibling" 
  US <- "unaffected sibling"
  for (event_freq in c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005)){
    #Adjust SNV freq
    CNV_SNV_table_proband_all <- CNV_SNV_table_proband_ori[CNV_SNV_table_proband_ori$event_freq <= event_freq, ]
    print(paste0("CNV_SNV_table_proband_all at ", event_freq*100, "% freq cutoff has ", nrow(CNV_SNV_table_proband_all), " rows"))
    CNV_SNV_table_unaff_sib_all <- CNV_SNV_table_unaff_sib_ori[CNV_SNV_table_unaff_sib_ori$event_freq <= event_freq, ]
    print(paste0("CNV_SNV_table_unaff_sib_all at ", event_freq*100, "% freq cutoff has ", nrow(CNV_SNV_table_unaff_sib_all), " rows"))

    ## group LOF variants in effect_priority
    CNV_SNV_table_proband_all$effect_priority[which(CNV_SNV_table_proband_all$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                          CNV_SNV_table_proband_all$typeseq_priority == "splicing")] <- "LoF" ###Bank
    CNV_SNV_table_unaff_sib_all$effect_priority[which(CNV_SNV_table_unaff_sib_all$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                          CNV_SNV_table_unaff_sib_all$typeseq_priority == "splicing")] <- "LoF" ###Bank                                      
    
    test.out <- data.frame()
    
    # for(freq in c(0.01, 0.001, 0.005, 0.0005, 0.0001)){ # event freq cut-offs
      for(pRecx in list(c(0,0.5), c(0.5,0.9),c(0.9, 1))){ # pRec cut-offs
        for(variant_type in c("damaging_missense", "LoF", "damaging_variant", "nonsynonymous", "synonymous")){
          for (GEgenes_only in c(TRUE, FALSE)) {
            if (GEgenes_only) {
              CNV_SNV_table_proband <- CNV_SNV_table_proband_all %>% subset(gene_symbol %in% GEgenes$V1)
              CNV_SNV_table_unaff_sib <- CNV_SNV_table_unaff_sib_all %>% subset(gene_symbol %in% GEgenes$V1)
            } else {
              CNV_SNV_table_proband <- CNV_SNV_table_proband_all
              CNV_SNV_table_unaff_sib <- CNV_SNV_table_unaff_sib_all
            }
            if (variant_type == "damaging_missense"){
              # transmitted variants in target in the proband and US tables
              target_proband.ix <- which(CNV_SNV_table_proband$gnomAD_pRec > pRecx[1] & CNV_SNV_table_proband$gnomAD_pRec <= pRecx[2] &
                                         CNV_SNV_table_proband$damaging_missense_count > 3)
              target_us.ix <- which(CNV_SNV_table_unaff_sib$gnomAD_pRec > pRecx[1] & CNV_SNV_table_unaff_sib$gnomAD_pRec <= pRecx[2] &
                                    CNV_SNV_table_unaff_sib$damaging_missense_count > 3)
            }
            if (variant_type == "LoF"){
              target_proband.ix <- which(CNV_SNV_table_proband$gnomAD_pRec > pRecx[1] & CNV_SNV_table_proband$gnomAD_pRec <= pRecx[2] &
                                         CNV_SNV_table_proband$effect_priority =="LoF")
              target_us.ix <- which(CNV_SNV_table_unaff_sib$gnomAD_pRec > pRecx[1] & CNV_SNV_table_unaff_sib$gnomAD_pRec <= pRecx[2] &
                                    CNV_SNV_table_unaff_sib$effect_priority =="LoF")
            }
            if (variant_type == "damaging_variant") {
              target_proband.ix <- which(CNV_SNV_table_proband$gnomAD_pRec > pRecx[1] & CNV_SNV_table_proband$gnomAD_pRec <= pRecx[2] &
                                         (CNV_SNV_table_proband$effect_priority =="LoF" | CNV_SNV_table_proband$damaging_missense_count > 3))
              target_us.ix <- which(CNV_SNV_table_unaff_sib$gnomAD_pRec > pRecx[1] & CNV_SNV_table_unaff_sib$gnomAD_pRec <= pRecx[2] &
                                    (CNV_SNV_table_unaff_sib$effect_priority =="LoF" | CNV_SNV_table_unaff_sib$damaging_missense_count > 3))
            }
            if (variant_type == "nonsynonymous"){
              target_proband.ix <- which(CNV_SNV_table_proband$gnomAD_pRec > pRecx[1] & CNV_SNV_table_proband$gnomAD_pRec <= pRecx[2] &
                                         CNV_SNV_table_proband$effect_priority =="nonsynonymous SNV")
              target_us.ix <- which(CNV_SNV_table_unaff_sib$gnomAD_pRec > pRecx[1] & CNV_SNV_table_unaff_sib$gnomAD_pRec <= pRecx[2] &
                                    CNV_SNV_table_unaff_sib$effect_priority =="nonsynonymous SNV")
            }
            if (variant_type == "synonymous"){
              target_proband.ix <- which(CNV_SNV_table_proband$gnomAD_pRec > pRecx[1] & CNV_SNV_table_proband$gnomAD_pRec <= pRecx[2] &
                                         CNV_SNV_table_proband$effect_priority =="synonymous SNV")
              target_us.ix <- which(CNV_SNV_table_unaff_sib$gnomAD_pRec > pRecx[1] & CNV_SNV_table_unaff_sib$gnomAD_pRec <= pRecx[2] &
                                    CNV_SNV_table_unaff_sib$effect_priority =="synonymous SNV")                     
            }
          
            # transmitted variants in proband in target
            target.proband <- nrow(unique(CNV_SNV_table_proband[intersect(target_proband.ix, which(CNV_SNV_table_proband$Relation == child)), 
                                                      c("#Sample", "#id")]))
            # transmitted variants in parent in target
            target.proband.parent <- sum(CNV_SNV_table_proband$Relation[target_proband.ix] == "parent")
            # nontransmitted variants in parent in target
            target.us.parent <- sum(CNV_SNV_table_unaff_sib$Relation[target_us.ix] == "parent")
            # nontransmitted variants in unaffected sibling in target
            target.us <- target.us.parent - nrow(unique(CNV_SNV_table_unaff_sib[intersect(target_us.ix, which(CNV_SNV_table_unaff_sib$Relation == US)), 
                                                      c("#Sample", "#id")]))
            target.child <- target.proband + target.us
            target.parent <- target.proband.parent + target.us.parent

            if (GEgenes_only) {
              target_proband.UID <- CNV_SNV_table_proband$UID[target_proband.ix]
              target_us.UID <- CNV_SNV_table_unaff_sib$UID[target_us.ix]

              target_proband_all.ix <- which(CNV_SNV_table_proband_all$UID %in% target_proband.UID)
              target_us_all.ix <- which(CNV_SNV_table_unaff_sib_all$UID %in% target_us.UID)

              # transmitted variants in parent not in target
              bg.proband <- sum(CNV_SNV_table_proband_all$Relation[-target_proband_all.ix] == child)
              # transmitted variants in parent not in target
              bg.proband.parent <- sum(CNV_SNV_table_proband_all$Relation[-target_proband_all.ix] == "parent")
              # nontransmitted variants in parent not in target
              bg.us.parent <- sum(CNV_SNV_table_unaff_sib_all$Relation[-target_us_all.ix] == "parent")
              # nontransmitted variants in US not in target
              bg.us <- bg.us.parent - sum(CNV_SNV_table_unaff_sib_all$Relation[-target_us_all.ix] == US)
              bg.child <- bg.proband + bg.us
              bg.parent <- bg.proband.parent + bg.us.parent
            } else {
              # transmitted variants in parent not in target
              bg.proband <- sum(CNV_SNV_table_proband$Relation[-target_proband.ix] == child)
              # transmitted variants in parent not in target
              bg.proband.parent <- sum(CNV_SNV_table_proband$Relation[-target_proband.ix] == "parent")
              # nontransmitted variants in parent not in target
              bg.us.parent <- sum(CNV_SNV_table_unaff_sib$Relation[-target_us.ix] == "parent")
              # nontransmitted variants in US not in target
              bg.us <- bg.us.parent - sum(CNV_SNV_table_unaff_sib$Relation[-target_us.ix] == US)
              bg.child <- bg.proband + bg.us
              bg.parent <- bg.proband.parent + bg.us.parent
            }

            df.test <- data.frame("target" = c(target.child, target.parent-target.child),
                                  "background" = c(bg.child, bg.parent-bg.child))
            test <- fisher.test(df.test, alternative = "greater")
            binom <- binom.test(c(df.test$target[1], df.test$target[2]))
            test.two.side <- fisher.test(df.test, alternative = "two.sided") 
            test.out <- rbind(test.out, 
                              data.frame("event_freq" = format(event_freq, scientific = F),
                                        "pRec" = pRecx[1],
                                        "variant_type" = variant_type,
                                        "gene_panel" = as.character(GEgenes_only),
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
      }
    #}
    table_name <- sprintf("%s_FisherTest_eventfreq%s.tsv", name, format(event_freq, scientific = F))
    write.table(test.out, table_name, sep="\t", row.names=F, quote=F, col.names=T)
  }
}

GEgenes <- fread("CHASE/recessive_genes/GE_Neurodevelopment_biallelic_Xavier_genes.txt", data.table = F, header = F)

CNV_SNV_table <- data.table::fread("path_to_output_from_step_3", data.table = F) 
nrow(CNV_SNV_table)
CNV_SNV_table <- CNV_SNV_table %>% subset(mosaic == FALSE) 
nrow(CNV_SNV_table)
CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$gnomAD_oe_lof_upper > 0.35 | is.na(CNV_SNV_table$gnomAD_oe_lof_upper)), ]
nrow(CNV_SNV_table)

all.meta <- fread("path_to_metadata_file", data.table = F)
names(all.meta) <- gsub(" ", ".", names(all.meta))
all.meta <- unique(all.meta)

meta_probandID <- all.meta$Sample.ID[all.meta$Relation %in% c('proband','affected sibling', 'child', 'sibling') &
                                       all.meta$Affection == 2]
meta_unaffSibID <- all.meta$Sample.ID[all.meta$Relation %in% c('child', 'other sibling', 'sibling', 'unaffected sibling') &
                                        all.meta$Affection == 1]

CNV_SNV_table_proband <- CNV_SNV_table %>% subset(sample %in% meta_probandID) 
CNV_SNV_table_unaff_sib <- CNV_SNV_table %>% subset(sample %in% meta_unaffSibID) 

nrow(CNV_SNV_table_proband)
nrow(CNV_SNV_table_unaff_sib)

Get_Fisher_Table(CNV_SNV_table_proband, CNV_SNV_table_unaff_sib, GEgenes, name = "parent_proband_US") # change name to desired output prefix
