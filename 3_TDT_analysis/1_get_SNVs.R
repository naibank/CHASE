library(data.table)

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
sample.path <- args[1]
snv_subset.outpath <- args[2]
print(args)

tmp <- fread(sample.path, data.table = FALSE)
tmp$"#Sample" <- as.character(tmp$"#Sample")
tmp <- tmp[which(tmp$high_quality == TRUE), 
             c("#Sample", "CHROM", "POS", "#id", "typeseq_priority", "effect_priority", "gene_symbol", "entrez_id", 
               "gene_type", "freq_max", "OZYG", "damaging_missense_count", "inheritance", "high_confidence_denovo",
               "gnomAD_pRec", "gnomAD_oe_lof_upper", "distance_spliceJunction")]
nrow(tmp)

if (!file.exists(snv_subset.outpath)) { 
    print("File doesn't exist, creating one")
    write.table(tmp, snv_subset.outpath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t") 
} else {
    print("File already exists, adding to it")
    write.table(tmp, snv_subset.outpath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", append = TRUE)
}







