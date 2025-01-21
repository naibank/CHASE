library(data.table)
library(dplyr)
options(echo=T)

# Functions ---------------------------------------------------------------

Get_CNV_SNV_Table <- function(cnv_sv, snvs, metadata) {
  table <- data.frame()
  for (cnv_row in 1:nrow(cnv_sv)) {
    message(sprintf("%s out of %s", cnv_row, nrow(cnv_sv)))
    cnv.tmp <- cnv_sv[cnv_row, names(cnv_sv)[!names(cnv_sv) %in% names(snvs)]]
    inh <- cnv.tmp$Inheritance
    Lof_cds_symbol <- strsplit(cnv.tmp$cds_symbol, "\\|")[[1]] #list of cds_symbols split at |
    
    probandID <- cnv.tmp$sample
    
    if (inh == "Paternal") {
      nonInhParentID <- metadata$Mother.ID[which(metadata$Sample.ID==probandID)]
      InhParentID <- metadata$Father.ID[which(metadata$Sample.ID==probandID)]
    } else {  #i.e. if inh == "Maternal"
      nonInhParentID <- metadata$Father.ID[which(metadata$Sample.ID==probandID)]
      InhParentID <- metadata$Mother.ID[which(metadata$Sample.ID==probandID)]
    }
    
    # SNVs that match gene
    snv <- snvs[snvs$gene_symbol %in% Lof_cds_symbol,]
    
    ## get SNVs of probands, transmitting and non-transmitting parent
    snv.inh.parent <- snv[which(snv$`#Sample` == InhParentID), ]
    snv.noninh.parent <- snv[which(snv$`#Sample` == nonInhParentID), ]
    snv.proband <- snv[which(snv$`#Sample` == probandID), ]
    
    ## get SNPs that found in both parents, as well as homozygous variants in non-inh parent. Then, remove them from the analysis
    snps.to.remove <- c(intersect(snv.inh.parent$"#id", snv.noninh.parent$"#id"),
                        snv.noninh.parent$"#id"[which(snv.noninh.parent$OZYG %in% c("alt-alt", "hom-alt"))])
    snv.noninh.parent <- snv.noninh.parent[which(!snv.noninh.parent$"#id" %in% snps.to.remove), ]
    
    ### retain only SNVs/Indels found in inh parent
    snv.proband <- snv.proband[which(snv.proband$"#id" %in% snv.noninh.parent$"#id"), ]
    snv.tmp <- rbind(snv.noninh.parent, snv.proband)
    
    ## check if snv.tmp is empty
    if (nrow(snv.tmp) == 0) {
      next #move to next cnv row
    }
    snv.tmp[, names(cnv.tmp)] <- cnv.tmp
    table <- bind_rows(table, snv.tmp)
  }
  return (table)
}


# Main --------------------------------------------------------------------

# snv.path <- "/Users/karahan/Desktop/CHASE/MSSNG_CG_all_snv.txt"
# meta.dir <- "/Users/karahan/Desktop/CHASE/metadata"
# cg.meta.path <- paste(meta.dir, "MSSNG_CG.metadata.tsv", sep = "/")
# ilmn.meta.path <- paste(meta.dir, "MSSNG_ILMN.metadata.tsv", sep = "/")
# ssc.meta.path <- paste(meta.dir, "SSC.metadata.tsv", sep = "/")
# spark1.meta.path <- paste(meta.dir, "SPARK_WGS_1.metadata.tsv", sep = "/")
# spark2.meta.path <- paste(meta.dir, "SPARK_WGS_2.metadata.tsv", sep = "/")
# spark3.meta.path <- paste(meta.dir, "SPARK_WGS_3.metadata.tsv", sep = "/")
# spark4.meta.path <- paste(meta.dir, "SPARK_WGS_4.metadata.tsv", sep = "/")
# cnv_sv.path <- "/Users/karahan/Desktop/CHASE/CNV+SV/MSSNG_CG.cnvs.svs.1per.cds.mosaic.tagged.tsv"

cg.snv.path <- "/hpf/largeprojects/tcagstor/scratch/kara.han/CHASE/data/SNV/MSSNG_CG_all_snv.txt"
ilmn.snv.path <- "/hpf/largeprojects/tcagstor/scratch/kara.han/CHASE/data/SNV/MSSNG_ILMN_all_snv.txt"
meta.dir <- "/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs"
cg.meta.path <- paste(meta.dir, "MSSNG_CG.metadata.tsv", sep = "/")
ilmn.meta.path <- paste(meta.dir, "MSSNG_ILMN.metadata.tsv", sep = "/")
cg.cnv_sv.path <- paste(meta.dir, "MSSNG_CG.cnvs.svs.10per.cds.mosaic.tagged.tsv", sep = "/")
ilmn.cnv_sv.path <- paste(meta.dir, "MSSNG_ILMN.cnvs.svs.10per.cds.mosaic.tagged.tsv", sep = "/")
outfile.path <- "/hpf/largeprojects/tcagstor/scratch/kara.han/CHASE/data/MSSNG_CNV_SNV_table.tsv"

cg.snvs <- fread(cg.snv.path, data.table = F)
ilmn.snvs <- fread(ilmn.snv.path, data.table = F)
snvs <- rbind(cg.snvs, ilmn.snvs)
names(snvs) <- c("#Sample", "CHROM", "POS", "#id", "typeseq_priority", "effect_priority", "gene_symbol", "entrez_id", 
                 "gene_type", "freq_max", "OZYG", "damaging_missense_count", "inheritance", "high_confidence_denovo",
                 "gnomAD_pRec", "gnomAD_oe_lof_upper", "distance_spliceJunction")
nrow(snvs)

# get trios only
cg.meta <- fread(cg.meta.path, data.table = F)
ilmn.meta <- fread(ilmn.meta.path, data.table = F)
all.meta <- rbind(cg.meta, ilmn.meta) 
names(all.meta) <- gsub(" ", ".", names(all.meta))
rm(cg.meta); rm(ilmn.meta)
nrow(all.meta)

trio.proband <- all.meta %>% subset(family_member == "trios") 
nrow(trio.proband)
sum(trio.proband$Mother.ID == "-")
sum(trio.proband$Father.ID == "-")
trio.proband <- trio.proband %>% subset(Mother.ID != "-" & Father.ID != "-") 
nrow(trio.proband)
sum(trio.proband$Mother.ID == "-")
sum(trio.proband$Father.ID == "-")
trio.all <- all.meta %>% subset(Sample.ID %in% c(trio.proband$Mother.ID, trio.proband$Father.ID, trio.proband$Sample.ID)) 
nrow(trio.all)
table(trio.all$Relation, trio.all$Affection)
parent02 <- trio.all %>% subset((Relation == "father" & Affection %in% c(0, 2)) | (Relation == "mother" & Affection %in% c(0, 2))) 
metadata <- trio.all %>% subset(! (Sample.ID %in% parent02$Sample.ID | Mother.ID %in% parent02$Sample.ID | Father.ID %in% parent02$Sample.ID)) 
nrow(parent02)
nrow(metadata)
table(metadata$Relation, metadata$Affection)

# get materally or paternally inherited CNVs only
cg.cnv_sv <- fread(cg.cnv_sv.path, data.table = F)
ilmn.cnv_sv <- fread(ilmn.cnv_sv.path, data.table = F)
cnv_sv <- rbind(cg.cnv_sv, ilmn.cnv_sv)
nrow(cnv_sv)
cnv_sv <- cnv_sv %>% subset(sample %in% metadata$Sample.ID)
cnv_sv <- cnv_sv %>% subset(Inheritance %in% c("Maternal", "Paternal"))
nrow(cnv_sv)

## get the list of genes and filter SNVs/Indels data based on the filtered list
cnv_sv.genes <- unique(strsplit(paste(cnv_sv$cds_symbol, collapse = "|"), "\\|")[[1]]) #list of cds_symbols split at |
snvs <- snvs %>% subset(gene_symbol %in% cnv_sv.genes)
nrow(snvs)

table <- Get_CNV_SNV_Table(cnv_sv, snvs, metadata)
nrow(table)

meta_parentsID <- metadata$Sample.ID[metadata$Relation %in% c('father','mother') & metadata$Affection == 1]
meta_probandID <- metadata$Sample.ID[metadata$Relation %in% c('proband','affected sibling', 'child', 'sibling') &
                                       metadata$Affection == 2]
meta_unaffSibID <- metadata$Sample.ID[metadata$Relation %in% c('child', 'other sibling', 'sibling', 'unaffected sibling') &
                                        metadata$Affection == 1]

## add a column to indicate the relation
table$Relation <- NA
table$Relation[which(table$`#Sample` %in% meta_parentsID)] <- "parent"
table$Relation[which(table$`#Sample` %in% meta_probandID)] <- "proband or affected sibling"
table$Relation[which(table$`#Sample` %in% meta_unaffSibID)] <- "unaffected sibling"
table <- table[which(!table$Relation == "NA"),]
table(table$Relation)

write.table(table, outfile.path, sep = "\t", row.names = F, quote = F, col.names = T)