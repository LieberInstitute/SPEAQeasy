#  This script is not intended to be run except from 'install_software.sh',
#  where the working directory is [SPEAQeasy path]/Software

print("Checking R packages...")

for (package in c("checkpoint", "here")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, repos = "http://cran.us.r-project.org")
    }
}

library("here")
library("checkpoint")

#  These are Bioconductor packages to install
packages <- c(
    "BiocParallel", "biomaRt", "Biostrings",
    "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Rnorvegicus.UCSC.rn6",
    "DelayedArray", "derfinder", "GenomicFeatures", "GenomicRanges",
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "rtracklayer",
    "SummarizedExperiment"
)

#  Install the Bioc packages and the GitHub package "jaffelab"
BiocManager::install(packages, update = FALSE)
remotes::install_github("LieberInstitute/jaffelab")

library("sessioninfo")
session_info()
