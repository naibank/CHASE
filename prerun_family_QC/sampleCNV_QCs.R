library(data.table)
library(GenomicRanges)
library(ggplot2)

#### perform sample QCs based on CNV count and SV count (3IQR+Q3 | 3SD+mean)
#### identify SVs/CNVs that are genomic disorders, impact ASD candidate genes, large >3MB and exclude those samples from the analysis
#### extract one proband and one unaffected sib per family
#### filter CNVs to CDS

# Steps:
# Remove cases with large CNVs > 3Mb, make sure to remove aneuploids
# Remove subjects with genomic syndromes, note the direction of CNV, deletion and duplication.
# Remove subjects with deletion impacting a gene on the genes in the list of clinically relevant CNVs.

datasets <- c(paste("SPARK_WGS_", 1:4, sep=""), "SSC", "MSSNG_ILMN", "MSSNG_CG")
cnv_files <- c(paste("CNVs.SPARK_WGS_", 1:4, ".freq_10percent.HQR.tsv", sep=""),
                "CNVs.SSC.freq_10percent.HQR.tsv",
                "CNVs.MSSNG_ILMN.freq_10percent.HQR.tsv",
                "CNVs.MSSNG_CG.freq_10percent.HQR.tsv")

sv_files <- c(paste("SVs.SPARK_WGS_", 1:4, ".freq_10percent.HQR.inh.tsv", sep=""),
              "SVs.SSC.freq_10percent.HQR.inh.tsv",
              "SVs.MSSNG_ILMN.freq_10percent.HQR.inh.tsv",
              "SVs.MSSNG_CG.freq_10percent.HQR.inh.tsv")

meta_files <- c(paste("SPARK_WGS_", 1:4, "_metadata.tsv", sep=""),
                "SSC_metadata.tsv",
                "MSSNG_metadata.tsv",
                "MSSNG_metadata.tsv")


# i <- 6
for(i in 1:7){
  cnvs <- fread(cnv_files[i], data.table = F)
  svs <- fread(sv_files[i], data.table = F)
  if("startAnn" %in% names(cnvs)){
    cnvs$STARTAnn <- cnvs$startAnn
    cnvs$ENDAnn <- cnvs$endAnn
  }
  
  if(datasets[i] == "MSSNG_CG"){
    svs$sv_type <- svs$Type
    svs$sv_type <- gsub("deletion", "DEL", svs$sv_type)
    svs[, c("chrAnn", "STARTAnn", "ENDAnn")] <- svs[, c("chrAnn", "startAnn", "endAnn")]
  }else{
    svs <- svs[which(svs$method == "DELLY|MANTA"), ] 
    svs[, c("chrAnn", "STARTAnn", "ENDAnn")] <- svs[, c("chrm", "start", "end")]
  }
  
  if(substr(svs$chrAnn[1], 1,3) != "chr"){
    svs$chrAnn <- paste0("chr", svs$chrAnn)
  }
  if(substr(cnvs$chrAnn[1], 1,3) != "chr"){
    cnvs$chrAnn <- paste0("chr", cnvs$chrAnn)
  }
  
  meta <- fread(meta_files[i], data.table = F)
  
  ### tag samples for removal from the analysis, removing big CNVs not SVs, as SVs not detected in CNVs are likely mosaic
  big_cnvs <- cnvs$sample[cnvs$SIZE > 3000000]
  big_cnvs <- data.frame("sample" = big_cnvs, "QC" = "big_CNVs")
  
  if(nrow(big_cnvs) > 0){
    cnvs <- cnvs[!cnvs$sample %in% big_cnvs$sample, ]
    svs <- svs[!svs$sample %in% big_cnvs$sample, ]
  }
  svs <- svs[svs$sv_type %in% c("DEL"), ]
  
  ### filter SVs that are smaller calls of a bigger call
  svs.g <- GRanges(svs$chrAnn, IRanges(svs$STARTAnn, svs$ENDAnn), "*")
  olap <- data.frame(findOverlaps(svs.g, svs.g))
  olap <- olap[olap$queryHits != olap$subjectHits, ]
  olap$width <- width(pintersect(svs.g[olap$queryHits], svs.g[olap$subjectHits]))
  olap$queryWidth <- olap$width/width(svs.g[olap$queryHits])
  olap$subjectWidth <- olap$width/width(svs.g[olap$subjectHits])
  olap <- olap[svs$sample[olap$queryHits] == svs$sample[olap$subjectHits] &
                 svs$sv_type[olap$queryHits] == svs$sv_type[olap$subjectHits], ]
  
  olap <- olap[olap$subjectWidth > 0.9, ]
  svs <- svs[-unique(olap$subjectHits), ]
  
  ### remove SVs already captured by CNVs, use SVs coordinates
  cnvs.g <- GRanges(cnvs$chrAnn, IRanges(cnvs$STARTAnn, cnvs$ENDAnn), "*")
  svs.g <- GRanges(svs$chrAnn, IRanges(svs$STARTAnn, svs$ENDAnn), "*")
  olap <- data.frame(findOverlaps(cnvs.g, svs.g))
  olap <- olap[cnvs$sample[olap$queryHits] == svs$sample[olap$subjectHits] &
                 cnvs$SVTYPE[olap$queryHits] == "DEL" & svs$sv_type[olap$subjectHits] == "DEL", ]
  
  olap$sample <- svs$sample[olap$subjectHits]
  olap$svs_start <- svs$startAnn[olap$subjectHits]
  olap$svs_end <- svs$endAnn[olap$subjectHits]
  
  svs_min_start <- aggregate(svs_start ~ queryHits, olap, min)
  svs_max_end <- aggregate(svs_end ~ queryHits, olap, max)
  svs_new_coordinate <- merge(svs_min_start, svs_max_end, by = "queryHits")
  
  olap <- merge(olap, svs_new_coordinate, by = "queryHits")
  svs <- svs[-unique(olap$subjectHits), ]
  
  olap <- olap[!duplicated(olap$queryHits), ]
  cnvs[unique(olap$queryHits), c("START", "END", "STARTAnn", "ENDAnn")] <- olap[, c("svs_start.y", "svs_end.y",
                                                                                    "svs_start.y", "svs_end.y")]
  ### qc based on cnvs
  del_count <- data.frame(table(cnvs$sample[cnvs$SVTYPE == "DEL"]))
  # dup_count <- data.frame(table(cnvs$sample[cnvs$SVTYPE == "DUP"]))
  
  cutoff_del <- min(3*IQR(del_count$Freq) + quantile(del_count$Freq, 0.75),
                3*sd(del_count$Freq) + mean(del_count$Freq))
  
  # cutoff_dup <- min(3*IQR(dup_count$Freq) + quantile(dup_count$Freq, 0.75),
  #                   3*sd(dup_count$Freq) + mean(dup_count$Freq))
  # 
  # cnv_qc <- unique(c(as.character(del_count$Var1[del_count$Freq > cutoff_del]),
  #               as.character(dup_count$Var1[dup_count$Freq > cutoff_dup])))
  cnv_qc <- del_count[del_count$Freq > cutoff_del, ] #2-1368-002
  cnv_qc$Freq <- "CNV_failed"
  names(cnv_qc) <- c("sample", "QC")
  
  ### qc based on svs
  del_count <- data.frame(table(svs$sample[svs$sv_type == "DEL"]))
  # dup_count <- data.frame(table(svs$sample[svs$sv_type == "DUP"]))
  
  cutoff_del <- min(3*IQR(del_count$Freq) + quantile(del_count$Freq, 0.75),
                    3*sd(del_count$Freq) + mean(del_count$Freq))
  
  # cutoff_dup <- min(3*IQR(dup_count$Freq) + quantile(dup_count$Freq, 0.75),
  #                   3*sd(dup_count$Freq) + mean(dup_count$Freq))
  # 
  # sv_qc <- unique(c(as.character(del_count$Var1[del_count$Freq > cutoff_del]),
  #                    as.character(dup_count$Var1[dup_count$Freq > cutoff_dup])))
  sv_qc <- del_count[del_count$Freq > cutoff_del, ]
  sv_qc$Freq <- "SV_failed"
  names(sv_qc) <- c("sample", "QC")
  
  # sv_count <- merge(del_count, dup_count, by = "Var1", all = T)
  # names(sv_count) <- c("Sample", "Del_count", "Dup_count")
  # sv_count$bhooma_qc <- FALSE
  # sv_count$bhooma_qc[sv_count$Sample %in% mssng_cnv_failed] <- TRUE
  # sv_count$rare_qc <- FALSE
  # sv_count$rare_qc[sv_count$Sample %in% sv_qc] <- TRUE
  # 
  # ggplot(sv_count, aes(x=bhooma_qc, y=Del_count)) + geom_boxplot(outlier.alpha = 0) + geom_point(aes(color = rare_qc), position = position_jitter(width = .2), alpha=.5) +
  #   scale_color_manual(values = c("grey", "red"))
  # 
  # ggplot(sv_count, aes(x=bhooma_qc, y=Dup_count)) + geom_boxplot(outlier.alpha = 0) + geom_point(aes(color = rare_qc), position = position_jitter(width = .2), alpha=.5) +
  #   scale_color_manual(values = c("grey", "red"))
  # 
  # sum(cnv_qc %in% mssng_cnv_failed)
  
  ### filter to cds impacted del
  cnvs <- cnvs[which(cnvs$cds_symbol != "" & cnvs$cds_symbol != "-" & !is.na(cnvs$cds_symbol)), ]
  svs <- svs[which(svs$cds_symbol != "" & svs$cds_symbol != "-" & !is.na(svs$cds_symbol)), ]
  svs$SVTYPE <- svs$sv_type
  
  
  if(i == 4){
    cnvs$del_freq_max  <- pmax(cnvs$cnvnIlmXPctFreq_50pctRecOvlp,
                           pmax(cnvs$cnvnIlm2PctFreq_50pctRecOvlp,
                                pmax(cnvs$erdsIlmXPctFreq_50pctRecOvlp,
                                     pmax(cnvs$erdsIlm2PctFreq_50pctRecOvlp, 
                                          pmax(cnvs$CGPctFreq_50pctRecOvlp,
                                               pmax(cnvs$otgCnvnPctFreq_50pctRecOvlp,
                                                    cnvs$otgErdsPctFreq_50pctRecOvlp)))))) 
  }else{
    cnvs$del_freq_max <- pmax(cnvs$cnvnIlmXParentalPercFreq_50percRecipOverlap,
                          pmax(cnvs$cnvnIlm2ParentalPercFreq_50percRecipOverlap,
                               pmax(cnvs$erdsIlmXParentalPercFreq_50percRecipOverlap,
                                    pmax(cnvs$erdsIlm2ParentalPercFreq_50percRecipOverlap, 
                                    pmax(cnvs$CGparentalPercFreq_50percRecipOverlap,
                                         pmax(cnvs$otgCnvnPercFreq_50percRecipOverlap,
                                              cnvs$otgErdsPercFreq_50percRecipOverlap))))))
  }
  
  
  if("CGparentalPercFreq_50percRecipOverlap" %in% names(svs)){
    svs$del_freq_max <- pmax(svs$CGparentalPercFreq_50percRecipOverlap,
                             pmax(svs$svMantaXPercFreq_90percRecipOverlap,
                                  pmax(svs$svManta2PercFreq_90percRecipOverlap,
                                       pmax(svs$svDellyXPercFreq_90percRecipOverlap,
                                            pmax(svs$svDelly2PercFreq_90percRecipOverlap,
                                                 pmax(svs$otgDellyPercFreq_90percRecipOverlap, 
                                                      svs$otgMantaPercFreq_90percRecipOverlap))))))
  }else if("svMantaXPercFreq_90percRecipOverlap" %in% names(svs)){
    svs$del_freq_max <- pmax(svs$CGparentalPercFreq_90percRecipOverlap,
                             pmax(svs$svMantaXPercFreq_90percRecipOverlap,
                                  pmax(svs$svManta2PercFreq_90percRecipOverlap,
                                       pmax(svs$svDellyXPercFreq_90percRecipOverlap,
                                            pmax(svs$svDelly2PercFreq_90percRecipOverlap,
                                                 pmax(svs$otgDellyPercFreq_90percRecipOverlap, 
                                                      svs$otgMantaPercFreq_90percRecipOverlap))))))
  }else{
    svs$del_freq_max <- pmax(svs$CGPctFreq_90pctRecOvlp,
                             pmax(svs$svMantaXPctFreq_90pctRecOvlp,
                                  pmax(svs$svManta2PctFreq_90pctRecOvlp,
                                       pmax(svs$svDellyXPctFreq_90pctRecOvlp,
                                            pmax(svs$svDelly2PctFreq_90pctRecOvlp,
                                                 pmax(svs$otgDellyPctFreq_90pctRecOvlp, 
                                                      svs$otgMantaPctFreq_90pctRecOvlp))))))
  }
  
  ### check 1kb or more svs
  # svs.1k <- svs[svs$sizeAnn > 1000, ]
  # cnvs.g <- GRanges(cnvs$CHROM, IRanges(cnvs$START, cnvs$END), "*")
  # svs.g <- GRanges(svs.1k$chrAnn, IRanges(svs.1k$start, svs.1k$end), "*")
  # 
  # olap <- data.frame(findOverlaps(cnvs.g, svs.g))
  # olap$sample <- cnvs$sample[olap$queryHits]
  # olap <- olap[olap$sample == svs.1k$sample[olap$subjectHits], ]
  # olap$width <- width(pintersect(svs.g[olap$subjectHits], cnvs.g[olap$queryHits]))/width(svs.g[olap$subjectHits])
  # olap <- olap[olap$width > 0.8, ]
  # 
  # svs.1k.notincnvs <- svs.1k[-olap$subjectHits, ]
  
  gsd <- read.delim("genomicSyndrome/Genomic.Syndromes.hg38.20200606.txt", stringsAsFactors = F)
  cnvs$pipeline <- "read-depth"
  svs$pipeline <- "paired-end"
  dt <- rbind(cnvs[, c("sample", "SVTYPE", "chrAnn", "STARTAnn", "ENDAnn", "exon_symbol", "variantTypeAnn", "cds_symbol", "pipeline", "Inheritance", "del_freq_max")],
              svs[, c("sample", "SVTYPE", "chrAnn", "STARTAnn", "ENDAnn", "exon_symbol", "variantTypeAnn", "cds_symbol", "pipeline", "Inheritance", "del_freq_max")])
  
  dt.g <- GRanges(dt$chrAnn, IRanges(dt$STARTAnn, dt$ENDAnn), "*")
  gsd.g <- GRanges(gsd$chr, IRanges(gsd$start, gsd$end), "*")
  olap <- data.frame(findOverlaps(dt.g, gsd.g))
  olap$overlap <- width(pintersect(dt.g[olap$queryHits], gsd.g[olap$subjectHits]))
  olap$overlap <- olap$overlap/width(gsd.g[olap$subjectHits])
  olap <- olap[olap$overlap > 0.9, ]
  olap$type <- ifelse(dt$SVTYPE[olap$queryHits] == "DEL", "Deletion", "Duplication")
  olap <- olap[olap$type == gsd$cnvFlag[olap$subjectHits], ]
  gd <- unique(dt$sample[olap$queryHits])
  
  crv <- readxl::read_excel("genomicSyndrome/MSSNG_SSC_pathogenic_SVs.xlsx")
  crv <- crv$Annotation
  crv <- crv[grep("del", crv)]
  crv <- gsub(" del", "", crv)
  
  inCRV <- function(genes, crv){
    return(sum(genes %in% crv) > 0) 
  }
  
  del <- dt[dt$SVTYPE == "DEL" & dt$cds_symbol != "" & !is.na(dt$cds_symbol), ]
  del.genes <- sapply(sapply(del$cds_symbol, strsplit, "\\|"), inCRV, crv)
  crv <- unique(del$sample[del.genes])
  
  explained <- unique(c(as.character(big_cnvs$sample), crv, gd))
  writeLines(explained, sprintf("%s.txt", datasets[i]))
  meta.f <- meta[!meta$`Sample ID` %in% c(explained, as.character(cnv_qc$sample), as.character(sv_qc$sample)), ]
  meta.f <- meta.f[meta.f$`Exclude due to Mendelian/gender problem?` == "no" & meta.f$`Exclude because re-sequenced?` == "no", ]
  
  proband.fam <- meta.f$`Family ID`[(!meta.f$Relation %in% c("father", "mother")) & meta.f$Affection == 2]
  unaff.fam <- meta.f$`Family ID`[!meta.f$Relation %in% c("father", "mother") & meta.f$Affection == 1]
  father.fam <- meta.f$`Family ID`[meta.f$Relation == "father"]
  mother.fam <- meta.f$`Family ID`[meta.f$Relation == "mother"]
  bothparent.fam <- intersect(father.fam, mother.fam)
  singleparent.fam <- setdiff(union(father.fam, mother.fam), bothparent.fam)
    
  probands <- meta.f[!meta.f$Relation %in% c("father", "mother") & meta.f$Affection == 2, ]
  probands <- probands[!duplicated(probands$`Family ID`), ]
  probands$family_member <- "singleton"
  probands$family_member[probands$`Family ID` %in% singleparent.fam] <- "dual"
  probands$family_member[probands$`Family ID` %in% bothparent.fam] <- "trios"
  
  
  unaff <- meta.f[!meta.f$Relation %in% c("father", "mother") & meta.f$Affection == 1, ]
  unaff <- unaff[!duplicated(unaff$`Family ID`), ]
  unaff$family_member <- "singleton"
  unaff$family_member[unaff$`Family ID` %in% singleparent.fam] <- "dual"
  unaff$family_member[unaff$`Family ID` %in% bothparent.fam] <- "trios"
  
  parent <- meta.f[meta.f$Relation %in% c("father", "mother"), ]
  parent <- parent[parent$`Family ID` %in% c(proband.fam, unaff.fam), ]
  parent$family_member <- ""
  
  meta <- rbind(probands, unaff, parent)
  if(datasets[i] == "MSSNG_ILMN")
    meta <- meta[meta$Platform != "Complete Genomics", ]
  if(datasets[i] == "MSSNG_CG")
    meta <- meta[meta$Platform == "Complete Genomics", ]
  
  meta$`Father ID`[meta$family_member == "singleton"] <- "-"
  meta$`Mother ID`[meta$family_member == "singleton"] <- "-"
  
  for(j in which(meta$family_member == "dual")){
    if(meta$`Mother ID`[j] %in% meta$Sample.ID){
      meta$`Father ID`[j] <- "-"
    }else if(meta$`Father ID`[j] %in% meta$Sample.ID){
      meta$`Mother ID`[j] <- "-"
    }
  }
  
  dt <- dt[dt$sample %in% meta$`Sample ID` & dt$SVTYPE == "DEL", ]
  
  write.table(dt, sprintf("%s.cnvs.svs.10per.cds.tsv", datasets[i]), sep="\t", row.names=F, quote=F, col.names=T)
  write.table(meta, sprintf("%s.metadata.tsv", datasets[i]), sep="\t", row.names=F, quote=F, col.names=T)
}
