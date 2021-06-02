#  Users at JHPCE will load an R module, with some system-wide packages
#  already installed, for components of the pipeline involving R scripts. This
#  script installs any missing R packages for the user.

packages <- c(
    "Biostrings", "biomaRt", "BiocParallel", "Biostrings",
    "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38",
    "BSgenome.Mmusculus.UCSC.mm10", "derfinder", "devtools",
    "GenomicFeatures", "GenomicRanges", "getopt",
    "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "plyr", "rafalib",
    "remotes", "rtracklayer", "SummarizedExperiment", "usethis"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

installed <- sapply(packages, requireNamespace, quietly = TRUE)

if (any(!installed)) {
    print(paste0("Missing the package(s) ", names(which(!installed))))

    BiocManager::install(packages[!installed], update = FALSE)
} else {
    print("All ordinary packages are already installed.")
}

#  Jaffelab is a special case since it is a package from GitHub
if (!("jaffelab" %in% installed.packages()[, 1])) {
    print("Missing the GitHub package 'jaffelab'.")
    library("remotes")
    remotes::install_github("LieberInstitute/jaffelab")
}

print("Done checking/installing packages.")

library("devtools")
devtools::session_info()
