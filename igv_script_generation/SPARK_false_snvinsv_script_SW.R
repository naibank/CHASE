setwd("/Users/shaniawu/Desktop/CHASE/SPARK_BA/IGV_SNV_check")

## Get SPARK metadata

SPARK_meta <- data.frame()

for (i in 1:3) {
  SPARK_metadata_path <- paste("../Data/SPARK_WGS_",i,"_metadata_relfixed.tsv",sep='')
  meta <- read.delim(SPARK_metadata_path, stringsAsFactors = F)
  meta$WGS_ver <- i # add col indicating SPARK WGS version number

  SPARK_meta <- rbind(SPARK_meta, meta)
}


## Get SNV table
IDs_and_pos <- fread('../ExtractSNVs/SPARK.parent.proband.SNVs.SNV1.CNV0.01.csv')

write_script <- function(text, path='igv_false_snvinsv_script.sh', 
                         append=T) {
  write(text, path, append=append)
}


for (i in 1:nrow(IDs_and_pos)) {
  write_script('new',)
  write_script('genome hg38')
  write_script('snapshotDirectory ./')
  
  row <- IDs_and_pos[i,]
  ID <- row$X.Sample
  chr <- row$CHROM
  pos <- row$POS
  
  # Get WGS version number for X.Sample
  ver <- SPARK_meta$WGS_ver[which(SPARK_meta$Sample.ID ==  row$X.Sample)]
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SPARK_WGS_%s/CRAM/%s.cram', ver,ID)
  write_script(sprintf('load %s',load_path))
  
  location <- sprintf('%s:%s-%s',chr,pos-100,pos+100)
  write_script(sprintf('goto %s',location))
  
  write_script('colorBy READ_STRAND')
  
  filename <- paste(ID,chr,pos,sep='_')
  write_script(sprintf('snapshot %s',paste(filename,'.png',sep='')))
}
