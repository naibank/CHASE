library(data.table)
load("/Users/shaniawu/SickKids CHASE/gene sets/gsMain_noLoFIn.Rdata")
CNV_SNV_table <- data.table::fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/CNV_SNV_table_eventfreq.tsv", data.table = F)

## group LOF variants in effect_priority
CNV_SNV_table$effect_priority[which(CNV_SNV_table$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") |
                                      CNV_SNV_table$typeseq_priority == "splicing")] <- "LoF" ###Bank


CNV_SNV_table <- CNV_SNV_table[which(CNV_SNV_table$Relation == "proband or affected sibling"), ] ###Bank - only use CH events in proband or affected sibling

all_genes <- unique(unlist(gsMain)) 
lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)
genes <- data.table::fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/hg38_refGene_20200708.transcript.txt", data.table = F)
genes <- unique(genes[, c("V5", "V6")])

table <- CNV_SNV_table[which(CNV_SNV_table$entrez_id %in% all_genes), ] 
  
freq <- 0.0001
pRecx <- 0.9

table <- table[which(table$event_freq < freq &
                        table$gene_symbol %in% lof$gene[lof$pRec > pRecx] &
                        table$effect_priority == "nonsynonymous SNV"), ]
#writeLines(unique(table$gene_symbol), "/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/MSSNG.nonsyn.0.01perc.0.9pRec.txt")

###load gsMain with LoF intorelant genes
load("/Users/shaniawu/SickKids CHASE/gene sets/gsMain_PGC_2021.RData")
all_genes <- unique(unlist(gsMain))
all_genes <- all_genes[which(!all_genes %in% table$entrez_id)] # remove genes in table

target.genes <- read.delim("/Users/shaniawu/SickKids CHASE/CNV_SNV_table by event_freq/data/MSSNG_genemania-genes50.txt", stringsAsFactors = F)
target.genes <- target.genes$Symbol
target.genes <- genes$V6[genes$V5 %in% target.genes] # get gene IDs
target.genes <- target.genes[target.genes %in% all_genes]

## make table
test.out <- data.frame()
  
for (gsi in 1:length(gsMain)){
    
  #Bank, get other genes for this frequency cut-off and variant type
  other_genes <- all_genes[which(!all_genes %in% target.genes)]
  
  target.CH <- sum(target.genes %in% gsMain[[gsi]])
  target.other <- sum(other_genes %in% gsMain[[gsi]])
  bg.CH <- sum(!target.genes %in% gsMain[[gsi]])
  bg.other <- sum(!other_genes %in% gsMain[[gsi]])
  df.test <- data.frame("target" = c(target.CH, target.other),
                                "background" = c(bg.CH, bg.other))
  test <- fisher.test(df.test, alternative = "greater")
  test.out <- rbind(test.out, 
                            data.frame("gene_set" = names(gsMain[gsi]), # gene set name
                                       "OR" = signif(test$estimate, digits = 3),
                                       "P" = signif(test$p.value, digits = 3),
                                       target.CH,
                                       target.other,
                                       bg.CH,
                                       bg.other))

}

test.out$BHFDR <- NA
test.out$BHFDR[grep("PhMm", test.out$gene_set)] <-
  p.adjust(test.out$P[grep("PhMm", test.out$gene_set)], method = "BH")
test.out$BHFDR[grep("Bspan|PhHs|Neurof|PSD|FMR1", test.out$gene_set)] <-
  p.adjust(test.out$P[grep("Bspan|PhHs|Neurof|PSD|FMR1", test.out$gene_set)], method = "BH")

table_name <- sprintf("MSSNG_CH_gene_enrichment_analysis_geneMANIA_50.tsv")
write.table(test.out, table_name, sep="\t", row.names=F, quote=F, col.names=T)

  

