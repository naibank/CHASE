#Required Packages
#BiocManager::install("GenomicRanges")
#BiocManager::install("Repitools")

library(GenomicRanges)
library(yaml)

##### ILMN DATA ##### (1111 CH events)

ILMNdata <- yaml::yaml.load_file("MSSNG_ILMN_CH_Data10P_Bank.yaml")
cnv <- data.frame(ILMNdata$AU2711303$CNVs)
cnv.g <- GRanges(cnv$CHROM, IRanges(cnv$START, cnv$END), "*")

genes <- data.table::fread("../../../ReferenceData/geneInfo2019/hg38_refGene_20200708.exon.txt", data.table = F)
genes <- genes[genes$V5 %in% strsplit(paste(cnv$gene_symbol, collapse = "|"), "\\|")[[1]], ]

genes.g <- GRanges(genes$V1, IRanges(genes$V2, genes$V3), "*")
genes.g <- reduce(genes.g)

olap <- data.frame(findOverlaps(cnv.g, genes.g))
olap$width <- width(pintersect(cnv.g[olap$queryHits], genes.g[olap$subjectHits]))
olap <- aggregate(width ~ queryHits, olap, sum)
cnv$exonicSize <- 0
cnv$exonicSize[olap$queryHits] <- olap$width

write.table(cnv, "cnv.with.exonicsize.tsv", sep="\t", row.names=F, quote=F, col.names=T)
