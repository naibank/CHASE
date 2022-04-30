setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")


IDs_and_pos <- fread('False SV check/MSSNG_SNV_falsecheck.csv')

write_script <- function(text, path='False SV check/igv_false_snvinsv_script.sh', 
                         append=T) {
  write(text, path, append=append)
}

for (i in 1:nrow(IDs_and_pos)) {
  write_script('new',)
  write_script('genome hg38')
  write_script('snapshotDirectory ./"IGV Results"/MSSNG_SNVinSV')
  
  row <- IDs_and_pos[i,]
  ID <- row$X.Sample
  chr <- row$CHROM
  pos <- row$POS
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/MSSNG/CRAM/%s.cram',ID)
  write_script(sprintf('load %s',load_path))
  
  location <- sprintf('%s:%s',chr,pos)
  write_script(sprintf('goto %s',location))
  
  write_script('colorBy READ_STRAND')
  
  filename <- paste(ID,chr,pos,sep='_')
  write_script(sprintf('snapshot %s',paste(filename,'.png',sep='')))
}

IDs_and_pos <- fread('False SV check/SSC_SNV_falsecheck.csv')

for (i in 1:nrow(IDs_and_pos)) {
  write_script('new',)
  write_script('genome hg38')
  write_script('snapshotDirectory ./"IGV Results"/SSC_SNVinSV')
  
  row <- IDs_and_pos[i,]
  ID <- row$X.Sample
  chr <- row$CHROM
  pos <- row$POS
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SSC/CRAM/%s.cram',ID)
  write_script(sprintf('load %s',load_path))
  
  location <- sprintf('%s:%s',chr,pos)
  write_script(sprintf('goto %s',location))
  
  write_script('colorBy READ_STRAND')
  
  filename <- paste(ID,chr,pos,sep='_')
  write_script(sprintf('snapshot %s',paste(filename,'.png',sep='')))
}