files <- list.files("../CH_Events/", full.names = T)

mssng <- files[grep("ILMN", files)]
ssc <- files[-grep("ILMN", files)]

dt.mssng <- data.frame()
for(file in mssng){
  suffix <- gsub(".tsv", "", strsplit(file, "_")[[1]][3])
  message(suffix)
  if(suffix %in% c("FATHER", "PATERNAL")){
    suffix <- "paternal"
  }else{
    suffix <- "maternal"
  }
  dt <- read.delim(file, stringsAsFactors = F)
  if(nrow(dt) > 0)
    dt$CNV_Inheritance <- suffix
  dt.mssng <- rbind(dt.mssng, dt)
}

dt.ssc <- data.frame()
for(file in ssc){
  suffix <- gsub(".tsv", "", strsplit(file, "_")[[1]][3])
  message(suffix)
  if(suffix %in% c("FATHER", "PATERNAL")){
    suffix <- "paternal"
  }else{
    suffix <- "maternal"
  }
  dt <- read.delim(file, stringsAsFactors = F)
  if(nrow(dt) > 0)
    dt$CNV_Inheritance <- suffix
  
  dt.ssc <- rbind(dt.ssc, dt)
}

mssng.meta <- read.delim("../../data/MSSNG_metadata.tsv", stringsAsFactors = F)
ssc.meta <- read.delim("../../data/SSC/SSC_metadata.tsv", stringsAsFactors = F)

mssng.meta <- mssng.meta[mssng.meta$Relation %in% c("affected sibling", "child", "proband"), c("Sample.ID", "Relation", "Mother.ID", "Father.ID")]
ssc.meta <- ssc.meta[ssc.meta$Relation %in% c("ohter sibling", "proband", "unaffected sibling"), c("Sample.ID", "Relation", "Mother.ID", "Father.ID")]

dt.mssng <- merge(dt.mssng, mssng.meta, by.x = "X.Sample", by.y = "Sample.ID", all = F)
dt.ssc <- merge(dt.ssc, ssc.meta, by.x = "X.Sample", by.y = "Sample.ID", all = F)

dt.mssng <- dt.mssng[dt.mssng$inheritance != "ambiguous" & dt.mssng$inheritance != dt.mssng$CNV_Inheritance, ]
dt.ssc <- dt.ssc[dt.ssc$inheritance != "ambiguous" & dt.ssc$inheritance != dt.ssc$CNV_Inheritance, ]

write.table(dt.mssng, "MSSNG.CHEvents.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(dt.ssc, "SSC.CHEvents.tsv", sep="\t", row.names=F, quote=F, col.names=T)
