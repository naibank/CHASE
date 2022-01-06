genes <- readLines("SSC.proband.nonsyn.0.01perc.0.9pRec.txt")
genes <- paste(genes, collapse="\t")

input <- readLines("sanitized_input.txt")
input[2] <- genes

writeLines(input, "SSC_probands_sanitized.txt")

genes <- readLines("SSC.unaffected_sibling.nonsyn.0.01perc.0.9pRec.txt")
genes <- paste(genes, collapse="\t")

input <- readLines("sanitized_input.txt")
input[2] <- genes

writeLines(input, "SSC_unaffectedsiblings_sanitized.txt")