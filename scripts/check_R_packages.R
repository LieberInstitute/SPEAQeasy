#  This script is not intended to be run except from 'install_software.sh',
#  where the working directory is [SPEAQeasy path]/Software

print("Checking R packages...")

lib_path <- paste0(getwd(), "/R-3.6.1/library")

if (!requireNamespace("checkpoint", quietly = TRUE)) {
      install.packages("checkpoint", repos = "http://cran.us.r-project.org", lib = lib_path)
  }

#  Automatically install ordinary packages as they existed when R 3.6.1 was
#  released
dir.create(paste0(getwd(), "/R-3.6.1/.checkpoint"))
checkpoint::checkpoint("2019-08-05",
    project = paste0(dirname(getwd()), "/scripts"),
    checkpointLocation = paste0(getwd(), "/R-3.6.1/")
)

#  These are Bioconductor packages to install
packages <- c(
    "BiocParallel", "biomaRt", "Biostrings",
    "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Rnorvegicus.UCSC.rn6",
    "DelayedArray", "derfinder", "GenomicFeatures", "GenomicRanges",
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "rtracklayer",
    "SummarizedExperiment"
)

#  Install the Bioc packages and the GitHub package "jaffelab"
BiocManager::install(packages, lib = lib_path, update = FALSE)
remotes::install_github("LieberInstitute/jaffelab", lib = lib_path)

library("sessioninfo")
session_info()
