#  Users at JHPCE will load an R module, with some system-wide packages
#  already installed, for components of the pipeline involving R scripts. This
#  script installs any missing R packages for the user.

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
    print(paste0("Missing the package(s) ", names(which(!installed))))

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
