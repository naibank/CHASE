library(data.table)

path <- "/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/exonic_splicing.ILMN"
files <- list.files(path, full.names = T)
total <- 2733
snvs <- data.frame()
for(i in 1696:total){
  message(sprintf("%s out of %s total variants = %s", i, total, nrow(snvs)))
  tmp <- fread(files[i], data.table = F)
  tmp$"#Sample" <- as.character(tmp$"#Sample")
  tmp <- tmp[which(tmp$high_quality == TRUE), 
             c("#Sample", "CHROM", "POS", "#id", "typeseq_priority", "effect_priority", "gene_symbol", "entrez_id", 
               "gene_type", "freq_max", "OZYG", "damaging_missense_count", "inheritance", "high_confidence_denovo",
               "gnomAD_pRec", "gnomAD_oe_lof_upper")]
  tmp$freq_max <- round(tmp$freq_max, digits = 3)
  snvs <- rbind(snvs, tmp)
  
  if(nrow(snvs) > 10000000){
    write.table(snvs, sprintf("./MSSNG.SNVs.freq._%s.tsv", i), sep="\t", row.names=F, quote=F, col.names=T)
    snvs <- data.frame()
  }
}
write.table(snvs, sprintf("./MSSNG.SNVs.freq._%s.tsv", total), sep="\t", row.names=F, quote=F, col.names=T)

