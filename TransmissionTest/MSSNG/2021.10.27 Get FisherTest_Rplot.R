library(data.table)
library(ggplot2)
setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")

metadata <- fread("MSSNG_metadata.tsv")

## find which IDs are parental and which are proband
meta_parentsID <- metadata$`Sample ID`[which(metadata$Relation == "mother" | metadata$Relation == "father")]
meta_probandID <- metadata$`Sample ID`[which(metadata$Relation == "proband")]
lof <- read.delim("/Users/shaniawu/SickKids CHASE/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)

res.out <- data.frame()
## separate table of probands and parents
for(freq in c(0.1, 0.05, 0.01, 0.001)){
  table <- fread("MSSNG_CNV_SNV_table.tsv")
  table <- table[which(!table$CHROM == "chrX")]
  table <- table[table$freq_max < freq, ]
  
  for(pRecx in c(0, 0.7, 0.8, 0.9)){
    table.prec <- table[table$gene_symbol %in% lof$gene[lof$pRec > pRecx], ]
    
    proband.ms <- sum(table.prec$"#Sample" %in% meta_probandID &
                        table.prec$effect_priority == "nonsynonymous SNV", na.rm = T)
    
    parent.ms <- sum(table.prec$"#Sample" %in% meta_parentsID &
                        table.prec$effect_priority == "nonsynonymous SNV", na.rm = T)
    
    
    proband.s <- sum(table.prec$"#Sample" %in% meta_probandID &
                        table.prec$effect_priority == "synonymous SNV", na.rm = T)
    
    parent.s <- sum(table.prec$"#Sample" %in% meta_parentsID &
                       table.prec$effect_priority == "synonymous SNV", na.rm = T)
    
    df.fet <- data.frame("ns" = c(proband.ms,parent.ms ), "s" = c(proband.s,parent.s))
    fet <- fisher.test(df.fet, alternative = "greater")
    
    res.out <- rbind(res.out,
                     data.frame(freq, pRecx, "OR" = fet$estimate,
                                "P" = fet$p.value))
  }
}


## plot -log(p-value) for each freq & pRec cut-offs 

res.out$neglog_P <- -log(res.out$P)
res.out$pRecx <- as.character(res.out$pRecx)
res.out$freq <- as.character(res.out$freq)

plot <- ggplot(data = res.out, aes(x=freq, y= neglog_P, fill=pRecx))
plot + geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "-log10(p-value)", x = "Frequency Cut-off") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "blue", size = 1)
  #facet_wrap(pRecx ~., nrow = 4)


