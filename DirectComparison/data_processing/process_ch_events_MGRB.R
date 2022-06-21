library(yaml)
library(survival)
library(GenomicRanges)
library(Repitools)
library(data.table)
library(dplyr)

setwd("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MGRB/data_processing/data/")

Determine_CNV_MAF <- function(CNVs) {
  CNV_MAF <- pmax(CNVs$CGparentalPercFreq_50percRecipOverlap,
                  pmax(CNVs$cnvnParentalPercFreq_50percRecipOverlap, 
                       pmax(CNVs$erdsParentalPercFreq_50percRecipOverlap,
                            na.rm=T), na.rm=T), na.rm=T)/100
  return (CNV_MAF)
}


Get_SNVs_With_Exonic_Size <- function(CNVs, SNVs, genes_data, CNV_freq_filter=1, 
                                      alt_base_filter_threshold=0.9) {
  # Determine CNV minor allele frequency
  total_exonic_size <- 0  
  return_SNV <- data.frame()
  if (length(CNVs) > 0) {
    CNVs$cnvFreq <- Determine_CNV_MAF(CNVs)
    CNVs <-  CNVs[which(CNVs$cnvFreq <= CNV_freq_filter), ]
    
    if (nrow(CNVs) > 0) {
      # # Rename startANN and endAnn cols in SVs to match CNVs:
      # if (!with_CNV) {
      #   names(CNVs)[names(CNVs) == 'startAnn'] <- 'STARTAnn'
      #   names(CNVs)[names(CNVs) == 'endAnn'] <- 'ENDAnn'
      # }
      
      CNVs_g <- GRanges(CNVs$chrAnn,
                        IRanges(CNVs$STARTAnn,
                                CNVs$ENDAnn), "*")
      genes <- genes_data[genes_data$V5 %in% strsplit(paste(CNVs$gene_symbol, collapse = "|"), "\\|")[[1]], ]
      
      if (nrow(genes) > 0) {
        genes_g <- GRanges(genes$V1, IRanges(genes$V2, genes$V3), "*")
        genes_g <- GenomicRanges::reduce(genes_g)
        
        olap <- data.frame(findOverlaps(CNVs_g, genes_g))
        
        if (nrow(olap) > 0) {
          olap$width <- width(pintersect(CNVs_g[olap$queryHits], genes_g[olap$subjectHits]))
          olap <- aggregate(width ~ queryHits, olap, sum)
          CNVs$exonicSize <- 0
          CNVs$exonicSize[olap$queryHits] <- olap$width
          total_exonic_size <- sum(CNVs$exonicSize)
        }
        if (nrow(SNVs) > 0) {
          SNVs_g <- GRanges(SNVs$X.CHROM, IRanges(SNVs$start, SNVs$end), "*")
          CNVs_g_redux <- GenomicRanges::reduce(CNVs_g)
          CNVs_g_redux_df <- data.frame(CNVs_g_redux)
          olap <- findOverlaps(SNVs_g, CNVs_g_redux)
          SNVs <- SNVs[unique(olap@from), ]
          SNVs$cnvchrAnn <- CNVs_g_redux_df$seqnames[olap@to]
          SNVs$cnvSTARTAnn <- CNVs_g_redux_df$start[olap@to]
          SNVs$cnvENDAnn <- CNVs_g_redux_df$end[olap@to]
          SNVs$cnvFreq <- CNVs$cnvFreq[olap@to]
          SNVs$CHfreq <- SNVs$cnvFreq * SNVs$freq_max
          SNVs$exonicSize <- CNVs$exonicSize[olap@to]#
          
          return_SNV <- rbind(return_SNV, SNVs)
        }
      }
    }
  }
  return_SNV <- subset(return_SNV, return_SNV$alt_fraction >= alt_base_filter_threshold)
  return(list(return_SNV, total_exonic_size))
}


Process_CNV_SNV_Data <- function(data, genes_data, 
                                 with_CNV = T, 
                                 CNV_freq_filter=1) {
  all_SNVs <- data.frame()
  all_exonic_sizes <- data.frame()
  all_CNVs <- data.frame()
  
  for(i in 1:length(data)){
    childCNVs <- data.frame(data[[i]]$CNVs) 
    childSNVs <- data.frame(data[[i]]$SNVs)
    
    if (nrow(childCNVs) > 0) {
      childData <- Get_SNVs_With_Exonic_Size(childCNVs, childSNVs, genes_data, 
                                             CNV_freq_filter=CNV_freq_filter)
      childSNVData <- childData[[1]]

      childExonicSize <- data.frame(Sample.ID=childCNVs$Sample.ID[1], exonicSize=childData[[2]])
      all_SNVs <- rbind(all_SNVs, childSNVData)
      all_exonic_sizes <- rbind(all_exonic_sizes, childExonicSize)
      all_CNVs <- dplyr::bind_rows(all_CNVs, childCNVs)
      all_CNVs$length <- (all_CNVs$END - all_CNVs$START) + 1 # add length column to CNVs
    }
    message(sprintf('%s of %s', i, length(data)))
  }
  return (list(all_SNVs, all_exonic_sizes, all_CNVs))
}


gene_exon_data <- data.table::fread("~/hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/gene_data/hg19_refGene_28_04_17.exon.txt", data.table = F)

CNV_data <- yaml::yaml.load_file("./MGRB_CH_Data_CNV1P_SNV.yaml")

CNV_SNVsExonicSizes_MGRB <- Process_CNV_SNV_Data(CNV_data,
                                                 gene_exon_data,
                                                 CNV_freq_filter=0.01)

write_yaml(CNV_SNVsExonicSizes_MGRB, "CNV_SNVsExonicSizes_MGRB.yaml")

