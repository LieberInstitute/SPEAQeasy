#  Install any missing R packages needed by SPEAQeasy

packages <- c(
    "BiocParallel", "biocthis", "biomaRt", "Biostrings",
    "checkpoint", "DelayedArray", "derfinder", "devtools", "GenomicFeatures",
    "GenomicRanges", "getopt", "here", "matrixStats", "org.Hs.eg.db",
    "org.Mm.eg.db", "org.Rn.eg.db", "plyr", "rafalib", "RColorBrewer",
    "rtracklayer", "sessioninfo", "styler", "SummarizedExperiment"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

installed <- sapply(packages, requireNamespace, quietly = TRUE)

if (any(!installed)) {
    print(paste0("Missing the package ", names(which(!installed))))

    BiocManager::install(packages[!installed], update = FALSE)
} else {
    print("All ordinary packages are already installed.")
}

#  Jaffelab is a special case since it is a package from GitHub
if (!requireNamespace("jaffelab", quietly = TRUE)) {
    print("Missing the GitHub package 'jaffelab'.")
    remotes::install_github("LieberInstitute/jaffelab")
}

print("Done checking/installing packages.")

devtools::session_info()
