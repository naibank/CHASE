library(data.table)

ls1 <- fread("/Users/shaniawu/SickKids CHASE/2021.11.01/data/CRVs/MSSNG+SSC.ASD135_LoF.tsv")
ls2 <- fread("/Users/shaniawu/SickKids CHASE/2021.11.01/data/CRVs/MSSNG+SSC.CNVs.tsv")
ls_exclude <- unique(rbind(ls1, ls2))

table <- fread("/Users/shaniawu/SickKids CHASE/CNV_SNV_table/CNV_SNV_table.tsv")

## add a column to indicate the relation
table$Relation <- NA
table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
table$Relation[which(table$`#Sample` %in% meta_probandID)] <- "proband or affected sibling"
table <- table[which(!table$Relation == "NA"),]
table <- table[which(!table$CHROM == "chrX")] #exclude variants on chrX


## get list of FamilyID & SampleID to exclude
metadata <- fread("/Users/shaniawu/Desktop/SickKids (F21 Coop)/CHASE (Bioinformatics Project)/2021.10.06/data/MSSNG_metadata.tsv")

familyID_exclude <- metadata$`Family ID`[which(metadata$`Sample ID` %in% ls_exclude$Sample)]
sampleID_exclude <- metadata$`Sample ID`[which(metadata$`Family ID` %in% familyID_exclude)]

## get new_table that excludes sampleID_exclude
new_table <- table[which(!table$`#Sample` %in% sampleID_exclude),]

write.table(new_table, "2021.11.01_CNV_SNV_table.tsv")


## plot table (old)
plot <- ggplot(data = table, aes(x = effect_priority, fill = Relation))
plot + geom_bar(stat = "count", position = position_dodge()) +
  geom_label(stat = "count", aes(label = ..count.., fill = Relation), position = position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot new_table
plot <- ggplot(data = new_table, aes(x = effect_priority, fill = Relation))
plot + geom_bar(stat = "count", position = position_dodge()) +
  geom_label(stat = "count", aes(label = ..count.., fill = Relation), position = position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

