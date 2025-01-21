
#### please check one parent sequence inherited CNVs
#### is it all CNVs with one parents, all only inherited CNVs that have "one_parent_sequenced" value

## STEP 1: Locate rare CNV deletions from ASD probands that have inherited it or P_denovo and envelop exons 
# Loading in CNV databases + changing first column name to Sample.ID 
# (necessary to keep identifying column name the same among all databases)
CNVfilter.procedure <- function(CNVfile.path, metatable){
  CNVs <- as_tibble(fread(CNVfile.path, data.table = FALSE))
  CNVs <- CNVs[!(CNVs$pipeline == "paired-end" & CNVs$mosaic), ]
  
  # To have uniform column name across all databases for Sample ID
  colnames(CNVs)[1] <- "Sample.ID"
  
  # Selecting only CNV deletions
  CNVs <- CNVs %>% dplyr::filter(variantTypeAnn == "del" | variantTypeAnn == "deletion")
  CNVs$size <- CNVs$ENDAnn - CNVs$STARTAnn + 1
  
  # qc cds del count
  cnv_count <- data.frame(table(CNVs$Sample.ID))
  cutoff <- min(3*IQR(cnv_count$Freq) + quantile(cnv_count$Freq, 0.75),
                3*sd(cnv_count$Freq) + mean(cnv_count$Freq))
  CNVs <- CNVs[which(!CNVs$Sample.ID %in% cnv_count$Var1[cnv_count$Freq > cutoff]), ]
  
  # denovo_count <- data.frame(table(CNVs$Sample.ID[which(CNVs$Inheritance == "P_denovo")]))
  
  #Filtering out CNVs with "NA", ambiguous and one parent sequenced inheritance
  parentCNVs <- CNVs  %>% dplyr::filter(is.na(Inheritance))
  #Process one parent CNVs; 
  for(one.i in which(CNVs$Inheritance == "One_parent_sequenced")){
    tmp.CNVs <- CNVs[one.i, ]
    tmp.sample <- metatable[metatable$Sample.ID == tmp.CNVs$Sample.ID, ]
    if(tmp.sample$Mother.ID != "-"){
      inh <- "Maternal"
      parent.ID <- tmp.sample$Mother.ID
    }else if(tmp.sample$Father.ID != "-"){
      inh <- "Paternal"
      parent.ID <- tmp.sample$Father.ID
    }
    
    tmp.parentCNVs <- parentCNVs[parentCNVs$Sample.ID == parent.ID, ]
    
    if(nrow(tmp.parentCNVs) > 0){
        tmp.CNVs <- GRanges(tmp.CNVs$chrAnn, IRanges(tmp.CNVs$STARTAnn, tmp.CNVs$ENDAnn), "*")
        tmp.parentCNVs <- GRanges(tmp.parentCNVs$chrAnn, IRanges(tmp.parentCNVs$STARTAnn, tmp.parentCNVs$ENDAnn), "*")
        olap <- data.frame(findOverlaps(tmp.CNVs, tmp.parentCNVs))
        if(nrow(olap) > 0){
          olap$overlap_width <- width(pintersect(
            tmp.CNVs[olap$queryHits], tmp.parentCNVs[olap$subjectHits]
          ))
          
          olap <- aggregate(overlap_width ~ queryHits, olap, sum)
          olap$overlap_width <- olap$overlap_width/width(tmp.CNVs[olap$queryHits])
          olap <- olap[olap$overlap_width > 0.5, ]
          
          if(nrow(olap) > 0){
            CNVs$Inheritance[one.i] <- inh
          }
        }
    }
  }
  
  #Selecting only Paternally or Maternally Inherited CNVs
  CNVs <- CNVs  %>% dplyr::filter(Inheritance %in% c("Paternal", "Maternal", "P_denovo", "Ambiguous", "Inherited_Ambiguous"))
  
  
  
  CNVs <- rbind(parentCNVs, CNVs)
  # Selecting CNVs which envelop at least one CDS region
  # CNVs <- CNVs %>% dplyr::filter(cds_symbol != ""& !is.na(cds_symbol)) ### the input is already cds
  
  # Selecting CNVs where the QC is "ok"
  # CNVs <- CNVs %>% dplyr::filter(Sample_QC == "ok") ### already did QCs in the preprocessing
  
  # Filtering out CNVs located on the sex chromosomes
  CNVs <- CNVs %>%  dplyr::filter(chrAnn != "chrX" & chrAnn != "chrY")
  # CNVs <- CNVs %>%  dplyr::filter(freq_max <= 1) ### the input is already 1 perc
  
  return(CNVs)
}

Get_CH_Events <- function(probands, CNVdf, SNVfolder){
  CH_Data <- list()
  # probands <- probands[probands$Sample.ID %in% CNVdf$Sample.ID, ]
  # write.csv(probands$Sample.ID, "TestSampleIDs_BE.csv")
  message(sprintf("Total %s probands", nrow(probands)))
  for(i in 1:nrow(probands)){
    
    # Get the filepath to the proband's relevant SNV file
    probandFilepath <- str_c(c(SNVfolder, "/", 
                               probands$Sample.ID[i], ".tsv.gz"), collapse = "") # Search proband SNV data
    
    proband <- probands$Sample.ID[i]
    
    # Empty vector for probands' mother
    mother <- vector("character", length = 1)
    
    # Get the filepath to mother SNV file
    mother <- probands %>% dplyr::filter(probands$Sample.ID == probands$Sample.ID[i]) %>% 
      pull(Mother.ID) # Search probands' mother SNV data
    motherFilepath <- str_c(c(SNVfolder, "/",
                              mother[1], ".tsv.gz"), collapse = "")
    
    # Empty vector for probands' father
    father <- vector("character", length = 1)
    
    # Get the filepath to father SNV file
    father <- probands %>% dplyr::filter(probands$Sample.ID == probands$Sample.ID[i]) %>%
      pull(Father.ID) #search probands' father SNV data
    fatherFilepath <- str_c(c(SNVfolder, "/",
                              father[1], ".tsv.gz"), collapse = "")
    
    # Skip probands with no parent ID information
    # if (probands$Sample.ID[i] == '-' | mother[1] == '-' | father[1] == '-') {  ###"""OLD
    #   #if ("-" %in% probands$Sample.ID | mother[1] == "-" | father[1] == "-") {
    #   
    #   next
    # }
    
    # Identify if proband has a CH event
    # METHOD 1: Identifies all SNVs that are within a CNV region of a proband (1) --> therefore, SNVs are on the allele without the deletion (CNV)
    # METHOD 2: Identifies all SNVs that are somewhat outside of CNV boundaries (2) --> however, located on the other allele + in a gene affected by the CNV
    
    ###"""WHERE ARE THE METHODS, CANT FIND SEPARATION BTWEEN 1 AND 2
    
    # Load CNV data from proband, mother and father
    probandCNVs <- CNVdf %>% dplyr::filter(Sample.ID == proband)
    motherCNVs <- CNVdf %>% dplyr::filter(Sample.ID == mother)
    fatherCNVs <- CNVdf %>% dplyr::filter(Sample.ID == father)
    
    CH_Data[[proband]]$CNVs <- probandCNVs
    CH_Data[[proband]]$Mother.CNVs <- motherCNVs
    CH_Data[[proband]]$Father.CNVs <- fatherCNVs

    # lens <- sapply(CH_Data[[proband]], nrow)
    
    # Get SNVs data
    for(member in c("proband", "Mother", "Father")){
      if(member == "proband"){
        cnv <- "CNVs"
        snvpath <- probandFilepath
        snvout <- "SNVs"
      }else if(member == "Mother"){
        cnv <- "Mother.CNVs"
        snvpath <- motherFilepath
        snvout <- "Mother.SNVs"
      }else{
        cnv <- "Father.CNVs"
        snvpath <- fatherFilepath
        snvout <- "Father.SNVs"
      }
      
      if(basename(snvpath) != "-.tsv.gz" & file.exists(snvpath)){
          snvs <- fread(snvpath, data.table = F)
          snvs <- snvs[snvs$high_quality == T, ]
          if(nrow(snvs) > 0 & nrow(CH_Data[[proband]][[cnv]]) > 0){
            
            cnvs.list <- unlist(sapply(CH_Data[[proband]][[cnv]]$exon_symbol, strsplit, "[|]"))
            duplicated.columns <- names(snvs)[duplicated(names(snvs))]
            
            snvs <- snvs[snvs$gene_symbol %in% cnvs.list, names(snvs)[!names(snvs) %in% duplicated.columns] ]
            
            if(nrow(snvs) > 0){
              snvs <- as_tibble(snvs) %>% mutate(LoF = ifelse(str_detect(effect_impact_str, "Stop_gain-High") == T | 
                                                                str_detect(effect_impact_str, "Frameshift-High") == T | 
                                                                str_detect(effect_impact_str, "Splice_site_High") == T |
                                                                str_detect(typeseq_priority, 'exonic;splicing') == T, T, F))
              CH_Data[[proband]][[snvout]] <- snvs  
            }else{
              CH_Data[[proband]][[snvout]] <- as_tibble(data.frame())
            }
            # snvs <- snvs[unique(olap@from), -c(222, 301)] #ILMN requires this
            
          }else{
            CH_Data[[proband]][[snvout]] <- as_tibble(data.frame())
          }
      }else{
        CH_Data[[proband]][[snvout]] <- as_tibble(data.frame())
      }
    }
    
    message(sprintf("%s, %s proband SNVs, %s father SNVs, %s mother SNVs", i, 
                    nrow(CH_Data[[proband]][["SNVs"]]), 
                    nrow(CH_Data[[proband]][["Father.SNVs"]]), 
                    nrow(CH_Data[[proband]][["Mother.SNVs"]])))
  }
  
  return(CH_Data)
}   

Get_SNVs <- function(file_path) {
  data <- yaml::yaml.load_file(file_path)
  
  childSNVs.in.CNVs <- data.frame()
  parentSNVs.in.parentCNVs <- data.frame()
  
  # childCount <- 0
  # fatherCount <- 0
  # motherCount <- 0
  for(i in 1:length(data)){
    
    childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
    childSNVs <- data.frame(data[[i]]$SNVs)
    
    fatherCNVs <- data.frame(data[[i]]$Father.CNVs)
    fatherSNVs <- data.frame(data[[i]]$Father.SNVs)
    
    motherCNVs <- data.frame(data[[i]]$Mother.CNVs)
    motherSNVs <- data.frame(data[[i]]$Mother.SNVs)
    
  
  
    # commonSNVs <- union(intersect(childSNVs$X.id, fatherSNVs$X.id),
    #                     intersect(childSNVs$X.id, motherSNVs$X.id))
    
    # if(nrow(childCNVs) > 0)
    #   childCNVs <- childCNVs[!childCNVs$mosaic & childCNVs$size > 1000, ]
    # if(nrow(fatherCNVs) > 0)
    #   fatherCNVs <- fatherCNVs[!fatherCNVs$mosaic & fatherCNVs$size > 1000, ]
    # if(nrow(motherCNVs) > 0)
    #   motherCNVs <- motherCNVs[!motherCNVs$mosaic & motherCNVs$size > 1000, ]
    
    # if(length(commonSNVs) > 0){
    #   childSNVs <- childSNVs[!childSNVs$X.id %in% commonSNVs, ]
    #   fatherSNVs <- fatherSNVs[!fatherSNVs$X.id %in% commonSNVs, ]
    #   motherSNVs <- motherSNVs[!motherSNVs$X.id %in% commonSNVs, ]
    # }
    #   childCount <- childCount + nrow(childSNVs)
    #   fatherCount <- fatherCount + nrow(fatherSNVs)
    #   motherCount <- motherCount + nrow(motherSNVs)
    # }
    
    if(nrow(childSNVs) > 0){
      childCNVs$cnvid <- paste(childCNVs$Sample.ID, childCNVs$chrAnn, childCNVs$STARTAnn, childCNVs$ENDAnn, childCNVs$Inheritance, sep="#")
      childCNVs$size <- childCNVs$ENDAnn - childCNVs$STARTAnn + 1
      childSNVs[, c("cnv_id", "cnv_chr", "cnv_start", "cnv_end", "size", "Inheritance", "pipeline", "del_freq_max")] <- NA
      for(snv_num in 1:nrow(childSNVs)){
        tmpcnv <- childCNVs[lapply(lapply(sapply(childCNVs$exon_symbol, strsplit, "\\|"), "==", childSNVs$gene_symbol[snv_num]), sum) > 0 , ]
        
        # tmpcnv <- childCNVs[which(childCNVs$chrAnn == childSNVs$CHROM[snv_num] &
        #                       childCNVs$STARTAnn < childSNVs$POS[snv_num] &
        #                       childCNVs$ENDAnn > childSNVs$POS[snv_num]), ]
        if(nrow(tmpcnv) > 0){
          childSNVs[snv_num, c("cnv_id", "cnv_chr", "cnv_start", "cnv_end", "size", "Inheritance", "pipeline", "del_freq_max")] <- 
            tmpcnv[, c("cnvid", "chrAnn", "STARTAnn", "ENDAnn", "size", "Inheritance", "pipeline", "del_freq_max")]
        }
      }  
      
      if(sum(is.na(childSNVs$cnv_id)) > 0)
        message("child with issue")
      childSNVs <- childSNVs[!is.na(childSNVs$cnv_id), ]
      
      if(nrow(childSNVs) > 0)
      if(nrow(childSNVs.in.CNVs) == 0){
        childSNVs.in.CNVs <- childSNVs
      }else{
        common.columns <- intersect(names(childSNVs.in.CNVs), names(childSNVs))
        childSNVs.in.CNVs <- rbind(childSNVs.in.CNVs[, common.columns], 
                                   childSNVs[, common.columns])
      }
    }
    
    if(nrow(fatherSNVs) > 0){
      fatherCNVs$cnvid <- paste(fatherCNVs$Sample.ID, fatherCNVs$chrAnn, fatherCNVs$STARTAnn, fatherCNVs$ENDAnn, fatherCNVs$Inheritance, sep="#")
      fatherCNVs$size <- fatherCNVs$ENDAnn - fatherCNVs$STARTAnn + 1
      fatherSNVs[, c("cnv_id", "cnv_chr", "cnv_start", "cnv_end", "size", "Inheritance", "pipeline", "del_freq_max")] <- NA
      for(snv_num in 1:nrow(fatherSNVs)){
        tmpcnv <- fatherCNVs[lapply(lapply(sapply(fatherCNVs$exon_symbol, strsplit, "\\|"), "==", fatherSNVs$gene_symbol[snv_num]), sum) > 0 , ]
        
        # tmpcnv <- fatherCNVs[which(fatherCNVs$chrAnn == fatherSNVs$CHROM[snv_num] &
        #                              fatherCNVs$STARTAnn < fatherSNVs$POS[snv_num] &
        #                              fatherCNVs$ENDAnn > fatherSNVs$POS[snv_num]), ]
        if(nrow(tmpcnv) > 0){
          fatherSNVs[snv_num, c("cnv_id", "cnv_chr", "cnv_start", "cnv_end", "size", "Inheritance", "pipeline", "del_freq_max")] <- tmpcnv[, c("cnvid", "chrAnn", "STARTAnn", "ENDAnn", "size", "Inheritance", "pipeline", "del_freq_max")]
        }
      }  
      
      if(sum(is.na(fatherSNVs$cnv_id)) > 0)
        message("father with issue")
      fatherSNVs <- fatherSNVs[!is.na(fatherSNVs$cnv_id), ]
      
      if(nrow(fatherSNVs) > 0)
      if(nrow(parentSNVs.in.parentCNVs) == 0){
        parentSNVs.in.parentCNVs <- fatherSNVs
      }else{
        common.columns <- intersect(names(parentSNVs.in.parentCNVs), names(fatherSNVs))
        parentSNVs.in.parentCNVs <- rbind(parentSNVs.in.parentCNVs[, common.columns], 
                                          fatherSNVs[, common.columns])
      }
    }
    
    if(nrow(motherSNVs) > 0){
      motherCNVs$cnvid <- paste(motherCNVs$Sample.ID, motherCNVs$chrAnn, motherCNVs$STARTAnn, motherCNVs$ENDAnn, motherCNVs$Inheritance, sep="#")
      motherCNVs$size <- motherCNVs$ENDAnn - motherCNVs$STARTAnn + 1
      motherSNVs[, c("cnv_id", "cnv_chr", "cnv_start", "cnv_end", "size", "Inheritance", "pipeline", "del_freq_max")] <- NA
      for(snv_num in 1:nrow(motherSNVs)){
        tmpcnv <- motherCNVs[lapply(lapply(sapply(motherCNVs$exon_symbol, strsplit, "\\|"), "==", motherSNVs$gene_symbol[snv_num]), sum) > 0 , ]
        
        # tmpcnv <- motherCNVs[which(motherCNVs$chrAnn == motherSNVs$CHROM[snv_num] &
        #                              motherCNVs$STARTAnn < motherSNVs$POS[snv_num] &
        #                              motherCNVs$ENDAnn > motherSNVs$POS)[snv_num], ]
        if(nrow(tmpcnv) > 0){
          motherSNVs[snv_num, c("cnv_id", "cnv_chr", "cnv_start", "cnv_end", "size", "Inheritance", "pipeline", "del_freq_max")] <- tmpcnv[, c("cnvid", "chrAnn", "STARTAnn", "ENDAnn", "size", "Inheritance", "pipeline", "del_freq_max")]
        }
      }  
      
      if(sum(is.na(motherSNVs$cnv_id)) > 0)
        message("mother with issue")
      motherSNVs <- motherSNVs[!is.na(motherSNVs$cnv_id), ]
      
      if(nrow(motherSNVs) > 0)
      if(nrow(parentSNVs.in.parentCNVs) == 0){
        parentSNVs.in.parentCNVs <- motherSNVs
      }else{
        common.columns <- intersect(names(parentSNVs.in.parentCNVs), names(motherSNVs))
        parentSNVs.in.parentCNVs <- rbind(parentSNVs.in.parentCNVs[, common.columns], 
                                          motherSNVs[, common.columns])
      }
    }
        
    message(sprintf("%s proband=%s, parent=%s (father=%s, mother=%s, proband=%s)", i, nrow(childSNVs.in.CNVs), nrow(parentSNVs.in.parentCNVs),
                   nrow(fatherSNVs), 
                   nrow(motherSNVs), 
                   nrow(childSNVs)))
    
  }
  
  return (rbind(childSNVs.in.CNVs, 
                parentSNVs.in.parentCNVs))
}


glmTest <- function(snvtable, metatable, gene_exon){
  if(nrow(snvtable[snvtable$X.Sample %in% metatable$Sample.ID, ]) > 0){
    snvtable <- snvtable[!duplicated(snvtable[, c("X.Sample", "X.id")]), ]
    
    # gene_exon.g <- GRanges(gene_exon$chr, IRanges(gene_exon$start, gene_exon$end), "*")
    # cnv.tmp <- unique(snvtable[, c("X.Sample", "cnv_chr", "cnv_start", "cnv_end")])
    # cnv.tmp.g <- GRanges(cnv.tmp$cnv_chr, IRanges(cnv.tmp$cnv_start, cnv.tmp$cnv_end), "*")
    # 
    # olap <- data.frame(findOverlaps(cnv.tmp.g, gene_exon.g))
    # olap$chr <- gene_exon$chr[olap$subjectHits]
    # olap$start <- gene_exon$start[olap$subjectHits]
    # olap$end <- gene_exon$end[olap$subjectHits]
    # olap$X.Sample <- cnv.tmp$X.Sample[olap$queryHits]
    # olap.g <- GRanges(paste(olap$X.Sample, olap$queryHits, sep="#"),
    #                   IRanges(olap$start, olap$end), "*")
    # olap.g <- data.frame(reduce(olap.g))
    # olap.g$seqnames <- as.character(olap.g$seqnames)
    # olap.g$X.Sample <- sapply(sapply(olap.g$seqnames, strsplit, "#"), "[", 1)
    # 
    # sizetable <- aggregate(width ~ X.Sample, olap.g, sum)
    # names(sizetable) <- c("X.Sample", "size")
    
    snvtable <- merge(snvtable, metatable[, c("Sample.ID", "Affection", "Sex", "Predicted.ancestry", "Dataset")], 
                      by.x = "X.Sample", by.y ="Sample.ID", all = F)
    names(gene_exon) <- c("X.Sample", "size")
    counttable <- data.frame(table(snvtable$"X.Sample"))
    counttable <- merge(counttable, gene_exon, by.y="X.Sample", by.x="Var1", all = T)
    
    metatable <- merge(metatable, counttable, by.x="Sample.ID", by.y = "Var1", all.x=T)
    metatable[is.na(metatable)] <- 0
    mod <- clogit(Affection ~ strata(Family.ID) + size + Freq, metatable)
    coeff <- mod$coefficients["Freq"]
    tmp <- summary(mod)
    p <- tmp$coefficients["Freq", "Pr(>|z|)"]
    
    kid <- sum(metatable$Affection == 2 & metatable$Freq > 0, na.rm = T)
    parent <- sum(metatable$Affection == 1 & metatable$Freq > 0, na.rm = T)
    
    return(data.table(kid, parent, coeff, p))
  }else{
    return(data.frame("kid" = NA, "parent" = NA, "coeff"=NA, "p"=NA))
  }
}
