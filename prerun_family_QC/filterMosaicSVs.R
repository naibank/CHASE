library(data.table)
library(GenomicRanges)

### remove mosaic deletions where SNVs found to be het instead of hom
args = commandArgs(trailingOnly=TRUE)

snv_files <- c("SPARK_WGS_1/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_2/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_3/variants/SNVs+indels/exonic+splicing/",
               "SPARK_WGS_4/variants/SNVs+indels/exonic+splicing/",
               "SSC/variants/SNVs+indels/exonic+splicing/",
               "MSSNG/ILMN/variants/SNVs+indels/exonic+splicing/",
               "MSSNG/CG/variants/SNVs+indels/exonic+splicing/")

datasets <- c(paste("SPARK_WGS_", 1:4, sep=""), "SSC", "MSSNG_ILMN", "MSSNG_CG")
i <- which(datasets == args[1])

out.cnvs <- data.frame()

cnvs <- read.delim(sprintf("%s.cnvs.svs.10per.cds.tsv", datasets[i]), stringsAsFactors = F)
uniq.samples <- unique(cnvs$sample)
uniq.samples <- uniq.samples[which(file.exists(paste(snv_files[i], uniq.samples, ".tsv.gz", sep="")))]

message(datasets[i])
for(sample in uniq.samples){
  message(sprintf("%s out of %s", which(uniq.samples == sample), length(uniq.samples)))
  tmp.cnvs <- cnvs[cnvs$sample == sample, ]
  snvs <- fread(paste(snv_files[i], sample, ".tsv.gz", sep=""), data.table = F)
  snvs <- snvs[which(snvs$high_quality), ]
  cnv.g <- GRanges(tmp.cnvs$chrAnn, IRanges(tmp.cnvs$STARTAnn, tmp.cnvs$ENDAnn), "*")
  snv.g <- GRanges(snvs$CHROM, IRanges(snvs$POS, snvs$POS), "*")
  
  olap <- data.frame(findOverlaps(snv.g, cnv.g))
  olap$zygosity <- snvs$OZYG[olap$queryHits]
  olap <- olap[olap$zygosity != "hom-alt", ]
  tmp.cnvs$mosaic <- F
  if(nrow(olap) > 0){
    tmp.cnvs$mosaic[unique(olap$subjectHits)] <- T
  }
  
  out.cnvs <- rbind(out.cnvs, tmp.cnvs)
}

write.table(out.cnvs, sprintf("%s.cnvs.svs.10per.cds.mosaic.tagged.tsv", datasets[i]), sep="\t", row.names=F, quote=F, col.names=T)
