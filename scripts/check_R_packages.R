#  This script is not intended to be run except from 'install_software.sh',
#  where the working directory is [SPEAQeasy path]/Software

print("Checking R packages...")

lib_path <- paste0(getwd(), "/R-4.1.2/library")

for (package in c("checkpoint", "here")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, repos = "http://cran.us.r-project.org", lib = lib_path)
    }
}

library("here")
library("checkpoint")

#  Automatically install ordinary packages as they existed a bit after R 4.1.2
#  was released
dir.create(here("Software", "R-4.1.2", ".checkpoint"))
checkpoint("2021-12-01",
    project_dir = here("scripts", "r_packages"),
    checkpoint_location = here("Software", "R-4.1.2")
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
