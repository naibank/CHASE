datasets <- c("MSSNG_CG", "MSSNG_ILMN", paste("SPARK_WGS_", 1:4, sep=""), "SSC")

if(!dir.exists("1_sh")){
  dir.create("1_sh")
}
sh <- c()
for(i in 1:length(datasets)){
  template <- readLines("run_template.sh")
  template[3] <- sprintf("#SBATCH --job-name %s", datasets[i])
  template[10] <- sprintf("Rscript 1_getCHEvent.R %s", datasets[i])
  
  writeLines(template, sprintf("1_sh/%s.sh", datasets[i]))
  sh <- c(sh, sprintf("sbatch 1_sh/%s.sh", datasets[i]))
}

writeLines(sh, "1_submit_jobs.sh")
