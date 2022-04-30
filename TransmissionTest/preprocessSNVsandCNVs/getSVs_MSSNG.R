library(data.table)

CNVfilter.procedure <- function(CNVfile.path){
  CNVs <- fread(CNVfile.path, stringsAsFactors = FALSE) 
  colnames(CNVs)[2] <- 'CHROM'
  colnames(CNVs)[3] <- 'START'
  colnames(CNVs)[4] <- 'END'
  colnames(CNVs)[132] <- 'Family.ID'
  colnames(CNVs)[133] <- 'Sample.ID'
  colnames(CNVs)[142] <- 'Mother.ID'
  colnames(CNVs)[143] <- 'Father.ID'
  # Selecting only CNV deletions
  CNVs <- CNVs[CNVs$sv_type == 'DEL', ]
  CNVs <- CNVs[CNVs$method == "MANTA" | CNVs$method == "DELLY|MANTA", ]
  # Filtering out CNVs with "NA", ambiguous and one parent sequenced inheritance
  # CNVs <- CNVs  %>% dplyr::filter(!is.na(Inheritance)) 
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance != "Ambiguous")
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance != "One_parent_sequenced") 
  
  # Selecting only Paternally or Maternally Inherited CNVs
  # CNVs <- CNVs  %>% dplyr::filter(Inheritance == "Paternal" | Inheritance == "Maternal")
  CNVs <- CNVs[CNVs$Inheritance %in% c('Paternal','Maternal'),]
  
  # Selecting CNVs which envelop at least one CDS region
  CNVs <- CNVs[CNVs$cds_symbol != "" & !is.na(CNVs$cds_symbol), ]
  
  # Selecting CNVs where the QC is "ok"
  CNVs <- CNVs[CNVs$Sample_QC == "PASS", ]
  
  # Filtering out CNVs located on the sex chromosomes
  CNVs <- CNVs[CNVs$CHROM != "chrX" & CNVs$CHROM != "chrY", ]
  
  # Filtering out CNVs with greater than MAF limit
  CNVs <- CNVs[CNVs$freq_max/100 <= 0.1,]
  
  return(CNVs)
}

tmp <- CNVfilter.procedure("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/MSSNG/SVs.ILMN/SVs.MSSNG_ILMN.freq_1percent.HQR.IntFreq.tsv")

write.table(tmp, "./MSSNG.SVs.freq.1perc.tsv", sep="\t", row.names=F, quote=F, col.names=T)

print('done')