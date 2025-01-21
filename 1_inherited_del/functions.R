
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
  
  denovo_count <- data.frame(table(CNVs$Sample.ID[which(CNVs$Inheritance == "P_denovo")]))
  
  #Filtering out CNVs with "NA", ambiguous and one parent sequenced inheritance
  parentCNVs <- CNVs  %>% dplyr::filter(is.na(Inheritance))
  #Process one parent CNVs; 
  # all CNVs from sample with one parent have an inheritance of "One_parent_sequenced"
  
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
  CNVs <- CNVs  %>% dplyr::filter(Inheritance %in% c("Paternal", "Maternal"))
  
  
  
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
  probands <- probands[probands$Sample.ID %in% CNVdf$Sample.ID, ]
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
    
    # get CNVs data per proband (mother and father CNVs as well)
    # if inheritance CNVs is not found in parent, we removed
    # filter out non transmitted CNVs
    # sometime parent might not have CNVs that transmitted to probands due to the fact that parents got filtered out during the QC steps 
    probandFilteredCNVs <- data.frame()
    for(inheritance in c("Maternal", "Paternal")){
      parentCNVs <- motherCNVs
      parent <- "Mother"
      if(inheritance == "Paternal"){
        parentCNVs <- fatherCNVs
        parent <- "Father"
      }
      
      if(sum(probandCNVs$Inheritance == inheritance, na.rm = T) > 0){
        inheritanceCNVs <- probandCNVs %>% dplyr::filter(Inheritance == inheritance)
        inheritanceCNVs.g <- GRanges(inheritanceCNVs$chrAnn, IRanges(inheritanceCNVs$STARTAnn, inheritanceCNVs$ENDAnn), "*")
        
        if(nrow(parentCNVs) > 0){
          parentCNVs.g <- GRanges(parentCNVs$chrAnn, IRanges(parentCNVs$STARTAnn, parentCNVs$ENDAnn), "*")
          olap <- data.frame(findOverlaps(inheritanceCNVs.g, parentCNVs.g))
          olap$width <- width(pintersect(inheritanceCNVs.g[olap$queryHits], parentCNVs.g[olap$subjectHits]))
          olap$maxwidth <- pmax(width(inheritanceCNVs.g[olap$queryHits]), width(parentCNVs.g[olap$subjectHits]))
          olap$rec_overlap <- olap$width/olap$maxwidth
          olap <- olap[olap$rec_overlap >= 0.5, ]
          if(nrow(olap) > 0){
            parentCNVs <- parentCNVs[unique(olap$subjectHits), ]
            probandFilteredCNVs <- rbind(probandFilteredCNVs, inheritanceCNVs[unique(olap$queryHits), ])
          }
        }
        
        CH_Data[[proband]][[sprintf("%s.CNVs", parent)]] <- parentCNVs
      }else{
        CH_Data[[proband]][[sprintf("%s.CNVs", parent)]] <- as_tibble(data.frame())
      }
    }
    
    CH_Data[[proband]]$CNVs <- probandFilteredCNVs
    
    # CH_Data[[proband]]$CNVs <- probandCNVs
    # CH_Data[[proband]]$Mother.CNVs <- motherCNVs
    # CH_Data[[proband]]$Father.CNVs <- fatherCNVs
    
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
  
  for(i in 1:length(data)){
    message(sprintf("%s proband=%s, parent=%s", i, nrow(childSNVs.in.CNVs), nrow(parentSNVs.in.parentCNVs)))
    childCNVs <- data.frame(data[[i]]$CNVs) # Could be proband or unaffected sibling
    childSNVs <- data.frame(data[[i]]$SNVs)
    
    fatherCNVs <- data.frame(data[[i]]$Father.CNVs)
    fatherSNVs <- data.frame(data[[i]]$Father.SNVs)
    
    motherCNVs <- data.frame(data[[i]]$Mother.CNVs)
    motherSNVs <- data.frame(data[[i]]$Mother.SNVs)
    
    #Bank started here
    if(nrow(childCNVs) > 0)
    for(cnv_num in 1:nrow(childCNVs)){
      if (!is.na(childCNVs$Inheritance[cnv_num])){
        if(childCNVs$Inheritance[cnv_num] == "Paternal"){
          parentCNVs <- fatherCNVs
          parentSNVs <- fatherSNVs
          
          ### remove parent SNVs found in proband
          snvtoremove <- union(parentSNVs$X.id[parentSNVs$X.id %in% childSNVs$X.id],
                               childSNVs$X.id[childSNVs$inheritance == "ambiguous"])
          parentSNVs <- parentSNVs[!parentSNVs$X.id %in% snvtoremove, ]
          
          tmpChildSNVs <- childSNVs[!is.na(childSNVs$inheritance) & childSNVs$inheritance != "paternal" & 
                                      childSNVs$inheritance != "ambiguous", ]
          tmpChildSNVs <- tmpChildSNVs[!tmpChildSNVs$X.id %in% snvtoremove, ]
        }else if (childCNVs$Inheritance[cnv_num] == "Maternal"){
          parentCNVs <- motherCNVs
          parentSNVs <- motherSNVs
          
          ### remove parent SNVs found in proband
          snvtoremove <- union(parentSNVs$X.id[parentSNVs$X.id %in% childSNVs$X.id],
                               childSNVs$X.id[childSNVs$inheritance == "ambiguous"])
          parentSNVs <- parentSNVs[!parentSNVs$X.id %in% snvtoremove, ]
          
          tmpChildSNVs <- childSNVs[!is.na(childSNVs$inheritance) & childSNVs$inheritance != "maternal" & 
                                      childSNVs$inheritance != "ambiguous", ]
          tmpChildSNVs <- tmpChildSNVs[!tmpChildSNVs$X.id %in% snvtoremove, ]
        }
      } else {
        break
      }
      
      childCNVs.g <- GRanges(childCNVs$chrAnn[cnv_num], 
                             IRanges(childCNVs$STARTAnn[cnv_num], childCNVs$ENDAnn[cnv_num]), "*")
      if(nrow(parentCNVs) > 0){
        parentCNVs.g <- GRanges(parentCNVs$chrAnn, IRanges(parentCNVs$STARTAnn, parentCNVs$ENDAnn), "*")
        
        olap <- data.frame(findOverlaps(childCNVs.g, parentCNVs.g))
        olap$width <- width(pintersect(childCNVs.g[olap$queryHits], parentCNVs.g[olap$subjectHits]))
        olap$maxwidth <- pmax(width(childCNVs.g[olap$queryHits]), width(parentCNVs.g[olap$subjectHits]))
        olap$rec_overlap <- olap$width/olap$maxwidth
        olap <- olap[olap$rec_overlap >= 0.5, ]
        
        if(nrow(olap) > 0){ ### if CNVs found parent, find CH events (SNVs that overlap CNVs in child and parent)
          #create cnv id
          cnvid <- paste(childCNVs$Sample.ID[cnv_num], childCNVs$chrAnn[cnv_num], childCNVs$STARTAnn[cnv_num], childCNVs$ENDAnn[cnv_num], childCNVs$Inheritance[cnv_num], sep="#")
          
          parentGenes <- strsplit(parentCNVs$exon_symbol[olap$subjectHits], ",")[[1]]
          tmpParentSNVs <- parentSNVs[parentSNVs$gene_symbol %in% parentGenes, ]
          
          ### get SNVs in parent
          # parentCNVs <- parentCNVs[olap$subjectHits, ]
          # parentCNVs.g <- parentCNVs.g[olap$subjectHits]
          if(nrow(tmpParentSNVs) > 0){
            #   parentSNVs.g <- GRanges(parentSNVs$CHROM, IRanges(parentSNVs$start, parentSNVs$end), "*")
            #   olap <- findOverlaps(parentSNVs.g, parentCNVs.g)
            #   
            #   if(length(olap) > 0){
            #     tmpParentSNVs <- parentSNVs[unique(olap@from), ]
            tmpParentSNVs$cnv_id <- cnvid
            tmpParentSNVs$cnv_chr <- childCNVs$chrAnn[cnv_num]
            tmpParentSNVs$cnv_start <- childCNVs$STARTAnn[cnv_num]
            tmpParentSNVs$cnv_end <- childCNVs$ENDAnn[cnv_num]
            tmpParentSNVs$size <- childCNVs$ENDAnn[cnv_num] - childCNVs$STARTAnn[cnv_num] + 1
            tmpParentSNVs$Inheritance <- childCNVs$Inheritance[cnv_num]
            tmpParentSNVs$del_freq_max <- childCNVs$del_freq_max[cnv_num]
            
            if(nrow(parentSNVs.in.parentCNVs) == 0){
              parentSNVs.in.parentCNVs <- tmpParentSNVs
              
            }else{
              common.columns <- intersect(names(parentSNVs.in.parentCNVs), names(tmpParentSNVs))
              parentSNVs.in.parentCNVs <- rbind(parentSNVs.in.parentCNVs[, common.columns], 
                                                tmpParentSNVs[, common.columns])
              
            }
          }
          # }
          ### get SNVs in kid
          childGenes <- strsplit(childCNVs$exon_symbol[cnv_num], ",")[[1]]
          tmpChildSNVs <- tmpChildSNVs[tmpChildSNVs$gene_symbol %in% childGenes, ]
          if(nrow(tmpChildSNVs) > 0){
            # tmpChildSNVs.g <- GRanges(tmpChildSNVs$CHROM, IRanges(tmpChildSNVs$start, tmpChildSNVs$end), "*")
            # olap <- findOverlaps(tmpChildSNVs.g, childCNVs.g)
            # 
            # if(length(olap) > 0){
            #   tmpChildSNVs <- tmpChildSNVs[unique(olap@from), ]
              tmpChildSNVs$cnv_id <- cnvid
              tmpChildSNVs$cnv_chr <- childCNVs$chrAnn[cnv_num]
              tmpChildSNVs$cnv_start <- childCNVs$STARTAnn[cnv_num]
              tmpChildSNVs$cnv_end <- childCNVs$ENDAnn[cnv_num]
              tmpChildSNVs$size <- childCNVs$ENDAnn[cnv_num] - childCNVs$STARTAnn[cnv_num] + 1
              tmpChildSNVs$Inheritance <- childCNVs$Inheritance[cnv_num]
              tmpChildSNVs$del_freq_max <- childCNVs$del_freq_max[cnv_num]
              
              if(nrow(childSNVs.in.CNVs) == 0){
                childSNVs.in.CNVs <- tmpChildSNVs
                
              }else{
                common.columns <- intersect(names(childSNVs.in.CNVs), names(tmpChildSNVs))
                childSNVs.in.CNVs <- rbind(childSNVs.in.CNVs[, common.columns], 
                                           tmpChildSNVs[, common.columns])
                
              }
            # }    
          }
          
        }else{
          ### do nothing
        }
      }
    }
  }
  
  return (rbind(childSNVs.in.CNVs, 
               parentSNVs.in.parentCNVs))
}


fisherTest <- function(snvtable, metatable){
  if(nrow(snvtable[snvtable$X.Sample %in% metatable$Sample.ID, ]) > 0){
    snvtable <- snvtable[!duplicated(snvtable[, c("X.Sample", "X.id")]), ]
    snvtable <- merge(snvtable, metatable[, c("Sample.ID", "Affection", "Sex", "Predicted.ancestry", "Dataset")], 
                      by.x = "X.Sample", by.y ="Sample.ID", all = F)

    counttable <- data.frame(table(snvtable$Affection, snvtable$cnv_id))
    
    counttable$Var2 <- as.character(counttable$Var2)
    counttable$Var1 <- as.numeric(counttable$Var1)
    counttable$Var1 <- ifelse(counttable$Var1 == 2, 1, 0)
    
    mod <- clogit(Var1 ~ strata(Var2) + Freq, counttable)
    coeff <- mod$coefficients
    tmp <- summary(mod)
    p <- tmp$coefficients[, "Pr(>|z|)"]
    
    kid <- sum(counttable$Freq[counttable$Var1 == 1])
    parent <- sum(counttable$Freq[counttable$Var1 == 0])
    
    return(data.table(kid, parent, coeff, p))
  }else{
    return(data.frame("kid"=NA, "parent"=NA, "coeff"=NA, "p"=NA))
  }
}

fisherTestProbandVUnaff <- function(probandSNV, probandtable, unaffSNV, unafftable){
  if(nrow(probandSNV[probandSNV$X.Sample %in% probandtable$Sample.ID, ]) > 0 & nrow(unaffSNV[unaffSNV$X.Sample %in% unafftable$Sample.ID, ]) > 0){
    
    probandSNV <- probandSNV[!duplicated(probandSNV[, c("X.Sample", "X.id")]), ]
    probandSNV <- merge(probandSNV, probandtable[, c("Sample.ID", "Affection", "Sex", "Predicted.ancestry", "Dataset")], 
                      by.x = "X.Sample", by.y ="Sample.ID", all = F)
    
    unaffSNV <- unaffSNV[!duplicated(unaffSNV[, c("X.Sample", "X.id")]), ]
    unaffSNV <- merge(unaffSNV, unafftable[, c("Sample.ID", "Affection", "Sex", "Predicted.ancestry", "Dataset")], 
                        by.x = "X.Sample", by.y ="Sample.ID", all = F)
    
    probandCount <- data.frame(table(probandSNV$Affection, probandSNV$cnv_id))
    unaffCount <- data.frame(table(unaffSNV$Affection, unaffSNV$cnv_id))
    
    continTable <- data.frame("X" = c(sum(probandCount$Freq[probandCount$Var1 == 2]),
                                      sum(unaffCount$Freq[unaffCount$Var1 == 2])),
                              "Y" = c(sum(probandCount$Freq[probandCount$Var1 == 1]),
                                      sum(unaffCount$Freq[unaffCount$Var1 == 1])))
    
    test <- fisher.test(continTable)
    coeff <- log(test$estimate)
    p <- test$p.value
    
    proband_ratio <- signif(continTable$X[1]/continTable$Y[1], digits = 2)
    unaff_ratio <- signif(continTable$X[2]/continTable$Y[2], digits = 2)
    
    return(data.table(proband_ratio, unaff_ratio, coeff, p))
  }else{
    return(data.frame("proband_ratio"=NA, "unaff_ratio"=NA, "coeff"=NA, "p"=NA))
  }
}
