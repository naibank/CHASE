################################################################################################################################################################################
# 
# MSSNG+SSC_process_CNV_SNV_table_SW.R
# purpose: combine MSSNG+SSC CNV_SNV_tables and add CH event frequency column 
# input: MSSNG_CNV_SNV_table.tsv, SSC_proband_CNV_SNV_table.tsv, SSC_unaffected_sibling_CNV_SNV_table.tsv
#         *database*.CNVs.freq.1perc.tsv, *database*.SVs.freq.1perc.tsv
#         /hpf/largeprojects/tcagstor/tcagstor_tmp/shania.wu/MSSNG+SSC/MSSNG+SSC_TDT/data
# output: MSSNG.SSC_CNV_SNV_table_eventfreq.tsv, SSC_unaffected_sibling_CNV_SNV_table_eventfreq.tsv
#             -> used as input for MSSNG+SSC_get_FisherTest_snvfreq_SW.R
# 
# notes: modified from MSSNG_process_CNV_SNV_table.R & SSC_process_CNV_SNV_table.R (Faraz)
#
##############################################################################################################################################################################

library(data.table)


##############################################################################################################################################################################
### MSSNG CNV_SNV_table
setwd("/hpf/largeprojects/tcagstor/scratch/kara.han/CHASE/data/10perc_CNV/CNV+SNV")

CNV_SNV_table <- fread("MSSNG_SSC_SPARK_CNV_SNV_table.tsv", data.table = F)
CNV_SNV_table$UID <- paste(CNV_SNV_table$`#Sample`, CNV_SNV_table$`#id`, CNV_SNV_table$sample, sep='.')
CNV_SNV_table <- dplyr::distinct(CNV_SNV_table, UID, .keep_all=T)
CNV_SNV_table$CNVID <- paste(CNV_SNV_table$sample, CNV_SNV_table$CHROM, CNV_SNV_table$START,
                             CNV_SNV_table$END)

cg.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/MSSNG_CG.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
ilmn.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/MSSNG_ILMN.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
ssc.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/SSC.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
spark1.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/SPARK_WGS_1.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
spark2.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/SPARK_WGS_2.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
spark3.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/SPARK_WGS_3.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
spark4.cnvs <- fread("/hpf/largeprojects/tcagstor/users/worrawat/CHASE/2024/prerun_1_sample_family_QCs/SPARK_WGS_4.cnvs.svs.10per.cds.mosaic.tagged.tsv", data.table=F)
cnvs <- rbind(cg.cnvs, ilmn.cnvs, ssc.cnvs, spark1.cnvs, spark2.cnvs, spark3.cnvs, spark4.cnvs)
colnames(cnvs)
colnames(cnvs)[11] <- 'CNV_freq'
cnvs$CNVID <- paste(cnvs$sample, cnvs$chrAnn, cnvs$STARTAnn, cnvs$ENDAnn)
## add new freq columns
CNV_SNV_table$CNV_freq <- NA
CNV_SNV_table$event_freq <- NA

# MSSNG.CNV_SNV_table <- merge(MSSNG.CNV_SNV_table, cnvs, by.x='sample', by.y='sample', all.x=F,all.y=F)
## add CNV freq column
for (row in 1:nrow(CNV_SNV_table)){
  message(row)
  CNV_SNV_table[row,]$CNV_freq <- cnvs$CNV_freq[cnvs$CNVID == CNV_SNV_table[row,]$CNVID][1]
}

## convert CNV percent to frequency
CNV_SNV_table$CNV_freq <- CNV_SNV_table$CNV_freq/100 

## add event_freq col
CNV_SNV_table$event_freq <- CNV_SNV_table$freq_max * CNV_SNV_table$CNV_freq ###Bank

write.table(CNV_SNV_table, "MSSNG_SSC_SPARK_CNV_SNV_table_eventfreq.tsv", sep="\t", row.names=F, quote=F, col.names=T)

print("done")


