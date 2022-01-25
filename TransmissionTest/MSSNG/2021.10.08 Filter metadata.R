library(data.table)

setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz")

metadata <- fread("./MSSNG_metadata.tsv", data.table = F)

Filter_Metadata <- function(file) {
  metadata <- fread(file)
  
  # List of samples belonging to families with a proband and both parents data:
  probands <- metadata$`Family ID`[which(metadata$Relation %in% c("proband", "affected sibling"))]
  mothers <- metadata$`Family ID`[which(metadata$Relation == "mother")]
  fathers <- metadata$`Family ID`[which(metadata$Relation == "father")]  
  families <- intersect(intersect(probands, mothers), fathers)  # a vector of Family IDs with a full family dataset
  
  samples_full_families <- metadata$`Sample ID`[which(metadata$`Family ID` %in% c(families))]
}

## Make sure that both parents have SNV data:

snv_1 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_1.tsv", data.table = F)
snv_2 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_2.tsv", data.table = F)
snv_3 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_3.tsv", data.table = F)
snv_4 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_4.tsv", data.table = F)
snv_5 <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG.SNVs.freq.10perc_5.tsv", data.table = F)

snv_1_sampleIDs <- snv_1$`#Sample`
snv_2_sampleIDs <- snv_2$`#Sample`
snv_3_sampleIDs <- snv_3$`#Sample`
snv_4_sampleIDs <- snv_4$`#Sample`
snv_5_sampleIDs <- snv_5$`#Sample`

all_snv_sampleIDs <- unique(union(union(union(union(snv_1_sampleIDs,snv_2_sampleIDs),
                                    snv_3_sampleIDs),snv_4_sampleIDs),snv_5_sampleIDs))

samples_full_families <- samples_full_families[which(samples_full_families %in% all_snv_sampleIDs)]

writeLines(samples_full_families, "2021.10.08_filtered_metadata.txt")

#### re-filter one more time
samples_full_families <- readLines("2021.10.08_filtered_metadata.txt")
metadata <- metadata[which(metadata$`Sample ID` %in% samples_full_families), ]

probands <- metadata$`Family ID`[which(metadata$Relation %in% c("proband", "affected sibling"))]
mothers <- metadata$`Family ID`[which(metadata$Relation == "mother")]
fathers <- metadata$`Family ID`[which(metadata$Relation == "father")]  
families <- intersect(intersect(probands, mothers), fathers)  # a vector of Family IDs with a full family dataset

metadata <- metadata[which(metadata$`Family ID` %in% families), ]

samples_full_families <- samples_full_families[which(samples_full_families %in% metadata$`Sample ID`)]
writeLines(samples_full_families, "2021.10.25_filtered_metadata_BE.txt")
