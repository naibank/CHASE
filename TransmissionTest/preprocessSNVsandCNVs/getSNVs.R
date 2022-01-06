library(data.table)

path <- "/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG/exonic_splicing.ILMN"
files <- list.files(path, full.names = T)

freq <- 0.1
count <- 1
snvs <- data.frame()
for(i in 1:length(files)){
  message(sprintf("%s out of %s total variants = %s", i, length(files), nrow(snvs)))
  tmp <- fread(files[i], data.table = F)
  tmp$"#Sample" <- as.character(tmp$"#Sample")
  tmp <- tmp[which(tmp$high_quality == TRUE & tmp$freq_max <= freq), 
             c("#Sample", "CHROM", "POS", "#id", "typeseq_priority", "effect_priority", "gene_symbol", "entrez_id", 
               "gene_type", "freq_max", "OZYG", "damaging_missense_count", "inheritance", "high_confidence_denovo")]
  tmp$freq_max <- round(tmp$freq_max, digits = 3)
  snvs <- rbind(snvs, tmp)
  
  if(i %% 2000 == 0){
    write.table(snvs, sprintf("../data/MSSNG.SNVs.freq.10perc_%s.tsv", count), sep="\t", row.names=F, quote=F, col.names=T)
    snvs <- data.frame()
    count = count + 1
  }
}
write.table(snvs, sprintf("../data/MSSNG.SNVs.freq.10perc_%s.tsv", count), sep="\t", row.names=F, quote=F, col.names=T)

