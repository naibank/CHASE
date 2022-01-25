

Get_CNVs_SNVs_From_YAML <- function(CH_events_file, CNV_freq_filter=0.1, SNV_freq_filter=0.1) {
  data <- yaml::yaml.load_file(CH_events_file)
  
  child_CNVs <- data.frame()
  child_SNVs <- data.frame()
  
  for (i in 1:length(data)) {
    curr_CNVs <- data.frame(data[[i]]$CNVs)
    if (nrow(curr_CNVs) > 0) {
      child_CNVs <- rbind(child_CNVs, curr_CNVs)
    }
    
    curr_SNVs <- data.frame(data[[i]]$SNVs)
    if (nrow(curr_SNVs) > 0) {
      names(curr_SNVs)[11] <- "SampleData"
      child_SNVs <- rbind(child_SNVs, curr_SNVs)
    }
    
    message(i)
  }
  
  child_CNVs <- child_CNVs[which(child_CNVs$freq_max < CNV_freq_filter), ]
  child_SNVs <- child_CNVs[which(child_SNVs$freq_max < SNV_freq_filter), ]
  
  return (list(child_CNVs, child_SNVs))
}

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz")
MSSNG_CNVs_SNVs <- Get_CNVs_SNVs_From_YAML("MSSNG_ILMN_CH_Data10P_Bank.yaml")
MSSNG_CNVs <- MSSNG_CNVs_SNVs[[1]]
MSSNG_SNVs <- MSSNG_CNVs_SNVs[[2]]

SSC_CNVs_SNVs <- Get_CNVs_SNVs_From_YAML("SSC_CH_Data10P_Bank.yaml")
SSC_CNVs <- SSC_CNVs_SNVs[[1]]
SSC_SNVs <- SSC_CNVs_SNVs[[2]]

write.table(MSSNG_CNVs, "MSSNG_CNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(MSSNG_SNVs, "MSSNG_SNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_CNVs, "SSC_CNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)
write.table(SSC_SNVs, "SSC_SNVs.tsv", sep = "\t", row.names = F, quote = F, col.names = T)


