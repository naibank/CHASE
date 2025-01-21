library(data.table)
library(survival)
library(GenomicRanges)
library(tidyverse)
library(stringr)
library(shiny)
library(scales)
library(tidyr)
source("functions.R")

ge <- readLines("../recessive_genes/GE_Neurodevelopment_biallelic_Xavier_genes.txt")

probandTable <- fread("SNV_in_CH_tables/Proband_parents.tsv", data.table = F)
unaffTable <- fread("SNV_in_CH_tables/Unaff_parents.tsv", data.table = F)

probandTable <- probandTable[which(probandTable$gnomAD_oe_lof_upper > 0.35 | is.na(probandTable$gnomAD_oe_lof_upper)), ]
unaffTable <- unaffTable[which(unaffTable$gnomAD_oe_lof_upper > 0.35 | is.na(unaffTable$gnomAD_oe_lof_upper)), ]

probandTable <- probandTable[probandTable$size > 20000, ]
unaffTable <- unaffTable[unaffTable$size > 20000, ]

gene_exon <- fread("../geneinfo/hg38_refGene_20200708.exon_headers.tsv", data.table = F)

metatable.all <- rbind(read.delim("../prerun_family_QCs/MSSNG_CG.metadata.tsv"),
                   read.delim("../prerun_family_QCs/MSSNG_ILMN.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_1.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_2.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_3.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SPARK_WGS_4.metadata.tsv"),
                   read.delim("../prerun_family_QCs/SSC.metadata.tsv"))

trios <- metatable.all[metatable.all$family_member == "trios" & metatable.all$Father.ID != "-" & metatable.all$Mother.ID != "-", ]


metatable <- rbind(trios, metatable.all[metatable.all$Sample.ID %in% c(trios$Father.ID, trios$Mother.ID), ])
affected_parent_fam <- metatable$Family.ID[metatable$Relation %in% c("father", "mother") & metatable$Affection %in% c(0, 2)]

metatable <- metatable[!metatable$Family.ID %in% 
                         affected_parent_fam, ]

### read CNVs
cnv_files <- list.files("../prerun_1_sample_family_QCs/", full.names = T,  pattern = "mosaic.tagged.tsv")
cnv_all <- data.frame()
for(f in cnv_files){
  cnv_all <- rbind(cnv_all, CNVfilter.procedure(f, metatable.all))
}

### constraint genes
cons.genes <- read.delim("../gene_info/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F,
                         header = T)

# gene_exon_size <- aggregate(size ~ genesymbol, gene_exon, sum)

cons.genes.list <- list()
cons.genes.list[["0"]] <- unique(gene_exon$genesymbol)
cons.genes.list[["0.5"]] <- unique(gene_exon$genesymbol[gene_exon$genesymbol %in% cons.genes$gene[which(cons.genes$pRec > 0.5)]])
cons.genes.list[["0.9"]] <- unique(gene_exon$genesymbol[gene_exon$genesymbol %in% cons.genes$gene[which(cons.genes$pRec > 0.9)]])
cons.genes.list[["ge"]] <- unique(gene_exon$genesymbol[gene_exon$genesymbol %in% ge])

### problem we had earlier, not all the exon were counted toward genome-wide correction, only those overlapping deletion were
#### get exon size table for correction
cnv_all <- data.frame(cnv_all)
cnv_size_table <- list()
for(i in 1:length(cons.genes.list)){
  exon.tmp <- gene_exon[gene_exon$genesymbol %in% cons.genes.list[[i]], ]
  exon.tmp.g <- GRanges(exon.tmp$chr, IRanges(exon.tmp$start, exon.tmp$end), "*")
  
  cnv_all.g <- GRanges(cnv_all$chrAnn, IRanges(cnv_all$STARTAnn, cnv_all$ENDAnn), "*")
  
  olap <- data.frame(findOverlaps(cnv_all.g, exon.tmp.g))
  olap$sample <- cnv_all$Sample.ID[olap$queryHits]
  olap$gene <- exon.tmp$genesymbol[olap$subjectHits]
  olap <- unique(olap[, c("sample", "gene")])
  
  olap <- merge(olap, exon.tmp, by.x = "gene", by.y = "genesymbol", all.x = T)
  olap.g <- GRanges(paste(olap$sample, olap$chr), IRanges(olap$start, olap$end), "*")
  olap.g <- GenomicRanges::reduce(olap.g)
  olap <- data.frame(olap.g)
  olap$sample <- sapply(sapply(as.character(olap$seqnames), strsplit, "\\ "), "[", 1)
  cnv_size_table[[names(cons.genes.list)[i]]] <- aggregate(width ~ sample, olap, sum)
}

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
    subtable <- probandTable#[probandTable$size > 10000, ]
    submeta <- metatable[!metatable$Relation %in% c("father", "mother") & metatable$Affection == 2, ]
    submeta <- rbind(submeta, metatable[metatable$Sample.ID %in% c(submeta$Father.ID, submeta$Mother.ID), ])
  } else{
    subtable <- unaffTable#unaffTable$size > 10000, ]
    submeta <- metatable[!metatable$Relation %in% c("father", "mother") & metatable$Affection == 1, ]
    submeta <- rbind(submeta, metatable[metatable$Sample.ID %in% c(submeta$Father.ID, submeta$Mother.ID), ])
    submeta$Affection[submeta$Relation %in% c("other", "other sibling", "child", "sibling", "unaffected sibling") &
                        submeta$Affection == 1] <- 2
  }
  
  for(var in c("synonymous SNV", "nonsynonymous SNV", "damaging missense", "LoF", "damaging variant")){
    if(var == "nonsynonymous SNV"){
      vartype = c(var, "damaging missense")
    }else if(var == "damaging variant"){
      vartype = c("LoF", "damaging missense")
    }else{
      vartype = var
    }
    
    for(freq in c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005)){
    for(pRec in list(c(0,0.5), c(0.5,0.9), c(0.9, 1))){
      
      gene_exon <- cnv_size_table[[as.character(pRec[1])]]
      
      res <- glmTest(subtable[which(subtable$effect_priority %in% vartype & 
                                      subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                         subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta, gene_exon)
           
      res.mssng <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & 
                                            subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                               subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta[submeta$Dataset == "MSSNG", ], gene_exon)
      
      res.ssc <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & 
                                          subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                               subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta[submeta$Dataset == "SSC", ], gene_exon)
      
      res.spark <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & 
                                            subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                               subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2]), ], submeta[grep("SPARK", submeta$Dataset), ], gene_exon)
      
      res.cnv <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & 
                                          subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                            subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2] &
                                          subtable$pipeline == "read-depth"), ], submeta, gene_exon)
      
      res.sv <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & 
                                         subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                          subtable$gnomAD_pRec > pRec[1] & subtable$gnomAD_pRec <= pRec[2] &
                                          subtable$pipeline == "paired-end"), ], submeta, gene_exon)
      
      names(res.mssng) <- paste("MSSNG", names(res.mssng), sep=".")
      names(res.ssc) <- paste("SSC", names(res.ssc), sep=".")
      names(res.spark) <- paste("SPARK", names(res.spark), sep=".")
      names(res.cnv) <- paste("RD", names(res.cnv), sep=".")
      names(res.sv) <- paste("PE", names(res.sv), sep=".")
      
      if(nrow(res) > 0){
        dt.out <- rbind(dt.out, 
                        data.frame(kid, var, freq, "pRec"=pRec[1], "gene_panel" = F, res, res.mssng[, 3:4], res.ssc[, 3:4], 
                                   res.spark[, 3:4], res.cnv[, 3:4], res.sv[, 3:4]))
      }             
  
      ###GE gene panel
      if(pRec[1] == 0){
        gene_exon <- cnv_size_table[["ge"]]
        
        res <- glmTest(subtable[which(subtable$effect_priority %in% vartype & subtable$gene_symbol %in% ge &
                                        subtable$del_freq_max/100 * subtable$freq_max < freq), ], submeta, gene_exon)
        
        res.mssng <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & subtable$gene_symbol %in% ge &
                                              subtable$del_freq_max/100 * subtable$freq_max < freq), ], submeta[submeta$Dataset == "MSSNG", ], gene_exon)
        
        res.ssc <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & subtable$gene_symbol %in% ge &
                                            subtable$del_freq_max/100 * subtable$freq_max < freq), ], submeta[submeta$Dataset == "SSC", ], gene_exon)
        
        res.spark <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & subtable$gene_symbol %in% ge &
                                              subtable$del_freq_max/100 * subtable$freq_max < freq), ], submeta[grep("SPARK", submeta$Dataset), ], gene_exon)
        
        res.cnv <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & subtable$gene_symbol %in% ge &
                                            subtable$del_freq_max/100 * subtable$freq_max < freq &
                                            subtable$pipeline == "read-depth"), ], submeta, gene_exon)
        
        res.sv <- glmTest(subtable[which(subtable$effect_priority  %in% vartype & subtable$gene_symbol %in% ge &
                                           subtable$del_freq_max/100 * subtable$freq_max < freq & 
                                           subtable$pipeline == "paired-end"), ], submeta, gene_exon)
        
        names(res.mssng) <- paste("MSSNG", names(res.mssng), sep=".")
        names(res.ssc) <- paste("SSC", names(res.ssc), sep=".")
        names(res.spark) <- paste("SPARK", names(res.spark), sep=".")
        names(res.cnv) <- paste("RD", names(res.cnv), sep=".")
        names(res.sv) <- paste("PE", names(res.sv), sep=".")
        
        if(nrow(res) > 0){
          dt.out <- rbind(dt.out, 
                          data.frame(kid, var, freq, "pRec"=pRec[1], "gene_panel" = T, 
                                     res, res.mssng[, 3:4], res.ssc[, 3:4], res.spark[, 3:4], res.cnv[, 3:4], res.sv[, 3:4]))
        }             
        
      }
      ################
    }
  }
}
}

dt.out[,1:8]
write.table(dt.out, "fisher.results.tsv", sep="\t", row.names=F, quote=F, col.names=T)
