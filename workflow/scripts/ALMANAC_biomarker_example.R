# Run: Rscript ALMANAC_biomarker_example.R /path/to/compat_v2_PSet.rds
# Keep ALMANAC_biomarker_helpers.R next to this example. Writes no files.
suppressPackageStartupMessages(library(PharmacoGx))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Supply the compatibility-v2 PSet RDS path")
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"ALMANAC_biomarker_helpers.R"))
p <- readRDS(args[1])
result <- almanac_pair_signature(p,c("Paclitaxel","Nilotinib"),
  features=c("ENSG00000000003","ENSG00000000005","ENSG00000000419"))
cat("Shared samples:",length(result$shared_samples),"\n")
cat("Excluded because RNA is absent:",paste(result$excluded_samples,collapse=", "),"\n")
print(result$signature)
