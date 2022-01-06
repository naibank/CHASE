genes <- readLines("MSSNG.nonsyn.0.01perc.0.9pRec.txt")
genes <- paste(genes, collapse="\t")

input <- readLines("sanitized_input.txt")
input[2] <- genes

writeLines(input, "MSSNG_sanitized.txt")