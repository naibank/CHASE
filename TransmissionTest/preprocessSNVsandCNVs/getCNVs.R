library(data.table)

tmp <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/alex.chan/MSSNG/CNVs.ILMN/CNVs.MSSNG_ILMN.freq_10percent.HQR.tsv", data.table = F)
tmp <- tmp[which(tmp$SVTYPE == "DEL" & tmp$Sample_QC == "ok" & tmp$Inheritance %in% c("Paternal", "Maternal") &
                 tmp$FILTER == "PASS"), 
           c("sample", "CHROM", "START", "END", "SVTYPE", "gene_symbol", "gene_egID", "exon_symbol", "exon_egID",
             "cds_symbol", "cds_egID", "dirtyRegion_percOverlap", "DGVpercFreq_subjects_coverageStudies_50percRecipOverlap",
             "cnvnIlmXParentalPercFreq_50percRecipOverlap", "cnvnIlm2ParentalPercFreq_50percRecipOverlap",
             "erdsIlmXParentalPercFreq_50percRecipOverlap", "erdsIlm2ParentalPercFreq_50percRecipOverlap",
             "sscErdsPercFreq_50percRecipOverlap", "sscCnvnPercFreq_50percRecipOverlap", "Inheritance")]

write.table(tmp, "../data/MSSNG.CNVs.freq.10perc.tsv", sep="\t", row.names=F, quote=F, col.names=T)
