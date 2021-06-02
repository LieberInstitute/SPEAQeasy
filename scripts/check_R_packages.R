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
    project = paste0(dirname(getwd()), "/scripts/"),
    checkpointLocation = paste0(getwd(), "/R-3.6.1/")
)

library("BiocManager")

#  These are Bioconductor packages to install
packages <- c(
    "Biostrings", "GenomicRanges", "GenomicFeatures", "org.Hs.eg.db",
    "biomaRt", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38",
    "org.Mm.eg.db", "BSgenome.Mmusculus.UCSC.mm10", "org.Rn.eg.db",
    "BSgenome.Rnorvegicus.UCSC.rn6", "bumphunter",
    "BiocParallel", "derfinder", "rafalib", "SummarizedExperiment",
    "plyr", "rtracklayer", "RColorBrewer", "LieberInstitute/jaffelab"
)

## Try to load the required packages and send a message for each package that fails to load
load_res <- sapply(packages, requireNamespace, quietly = TRUE)
if (any(!load_res)) {
    print(paste0("check_R_packages.R: missing the package: ", names(which(!load_res))))

    #  Install missing packages: suppressing updates to avoid "path not writeable" warnings
    ## when trying to update packages above the user's file permissions
    present <- installed.packages(lib.loc = lib_path)
    BiocManager::install(packages[!packages %in% rownames(present)],
        update = TRUE,
        lib = lib_path
    )
} else {
    print("check_R_packages.R: Everything OK. No R packages need to be installed")
}

library("devtools")
devtools::session_info()
