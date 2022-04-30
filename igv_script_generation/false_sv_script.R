setwd("/Users/faraz/Documents/Work/Workterm 5 - SickKids/CHASE Project/Faraz/Data")


IDs_and_pos <- fread('False SV check/SSC_SV_1kto2.5k.csv')

write_script <- function(text, path='False SV check/igv_falsesv_ssc_script.sh', 
                         append=T) {
  write(text, path, append=append)
}

for (i in 1:nrow(IDs_and_pos)) {
  write_script('new',)
  write_script('genome hg38')
  write_script('snapshotDirectory ./')
  
  row <- IDs_and_pos[i,]
  ID <- row$sample
  start <- row$START
  end <- row$END
  chr <- row$chrAnn
  
  start <- start-1000
  end <- end+1000
  
  load_path <- sprintf('https://bounce2.tcag.ca/tcag/mirrored/SSC/CRAM/%s.cram',ID)
  write_script(sprintf('load %s',load_path))
  
  location <- sprintf('%s:%s-%s',chr,start,end)
  write_script(sprintf('goto %s',location))
  
  write_script('viewaspairs')
  
  write_script('snapshot')
}