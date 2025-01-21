datasets <- c(paste("SPARK_WGS_", 1:4, sep=""), "SSC", "MSSNG_ILMN", "MSSNG_CG")

if(!dir.exists("2_sh")){
  dir.create("2_sh")
}
sh <- c()
for(i in 1:length(datasets)){
  template <- readLines("run_template.sh")
  template[3] <- sprintf("#SBATCH --job-name %s", datasets[i])
  template[10] <- sprintf("Rscript filterMosaicSVs.R %s", datasets[i])
  
  writeLines(template, sprintf("2_sh/%s.sh", datasets[i]))
  sh <- c(sh, sprintf("sbatch 2_sh/%s.sh", datasets[i]))
}

writeLines(sh, "submit_jobs.sh")
