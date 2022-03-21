#  This script is not intended to be run except from 'install_software.sh',
#  where the working directory is [SPEAQeasy path]/Software

print("Checking R packages...")

# base R packages:
plist <- c("tidyverse", "data.table", "devtools", "sqldf", "remotes", "stringr", 
   "stringi", "usethis", "here", "reshape2")
for (p in plist){
 if(! p %in% installed.packages()) install.packages(p, dependencies = TRUE, type = "both")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

#  These are Bioconductor packages to install
packages <- c(
    "BiocParallel", "biomaRt", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Rnorvegicus.UCSC.rn6",
    "Matrix", "DelayedArray", "derfinder", "GenomicFeatures", "GenomicRanges",
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "rtracklayer", "limma",
    "SummarizedExperiment"
)

#  Install the Bioc packages and the GitHub package "jaffelab"
BiocManager::install(packages, update = FALSE)
if(! 'jaffelab' %in% installed.packages())
     devtools::install_github('LieberInstitute/jaffelab')
