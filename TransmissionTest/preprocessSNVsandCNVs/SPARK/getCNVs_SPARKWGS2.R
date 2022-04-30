library(data.table)

tmp <- fread("/hpf/largeprojects/tcagstor/tcagstor_tmp/farazali/SPARK_WGS_2/CNVs/CNVs.SPARK_WGS_2.freq_1percent.HQR.tsv", data.table = F)
tmp <- tmp[which(tmp$SVTYPE == "DEL" & tmp$Sample_QC == "pass" & tmp$Inheritance %in% c("Paternal", "Maternal") &
                   tmp$FILTER == "PASS"), 
           c("sample", "CHROM", "START", "END", "SVTYPE", "gene_symbol", "gene_egID", "exon_symbol", "exon_egID",
             "cds_symbol", "cds_egID", "dirtyRegion_percOverlap", "DGVpercFreq_subjects_coverageStudies_50percRecipOverlap",
             "CGparentalPercFreq_50percRecipOverlap", "otgErdsPercFreq_50percRecipOverlap", "otgCnvnPercFreq_50percRecipOverlap",
             "cnvnIlmXParentalPercFreq_50percRecipOverlap", "cnvnIlm2ParentalPercFreq_50percRecipOverlap",
             "erdsIlmXParentalPercFreq_50percRecipOverlap", "erdsIlm2ParentalPercFreq_50percRecipOverlap",
             "sscErdsPercFreq_50percRecipOverlap", "sscCnvnPercFreq_50percRecipOverlap", "Inheritance")]

write.table(tmp, "./SPARKWGS2.CNVs.freq.1perc.tsv", sep="\t", row.names=F, quote=F, col.names=T)
