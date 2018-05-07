## Start an outmess object to store messages
outmess <- "Log: checking R packages\n=============="

## Check if devtools package is installed
## If it is not installed, install it
if (!requireNamespace("devtools", quietly = TRUE)) {

  outmess <- c(outmess, "installing devtools package")

  source("https://bioconductor.org/biocLite.R")
  ## Suppresing updates to avoid "path not writeable" warnings
  biocLite('devtools', suppressUpdates = T)
}

# Load devtools package
library('devtools')

## Check if jaffelab package is installed
## If it is not installed, install it
if (!requireNamespace("jaffelab", quietly = TRUE)) {

  outmess <- c(outmess,"installing jaffelab package")

  install_github('LieberInstitute/jaffelab')
}

# Set the vector with package names required by the pipeline
# devtools and jaffelab were removed from this install block since they are installed in the previous block
packages <- c('Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db',
              'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
              'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db',
              'BSgenome.Rnorvegicus.UCSC.rn6', 'derfinder', 'bumphunter',
              'getopt', 'BiocParallel', 'rafalib', 'SummarizedExperiment',
              'plyr', 'rtracklayer', 'RColorBrewer')

## Try to load them
load_res <- sapply(packages, requireNamespace, quietly = TRUE)

## If any of the packages could not be loaded, send a message
if(any(!load_res)) {
  outmess <- c(outmess,paste0("check_R_packages.R: missing the package: ",names(which(load_res == F))))
}

## Store all the packages that are not installed
not_installed_packages <- names(which(load_res == F))

## Check if there are packages to be installed and procced accordingly
if( length(not_installed_packages) == 0) {

  outmess <- c(outmess,"check_R_packages.R: Everything OK. No R packages need to be installed")

} else {

  outmess <- c(outmess,"check_R_packages.R: Attempting to install missing packages")

  ## Install missing packages
  source('https://bioconductor.org/biocLite.R')

  ## Suppresing updates to avoid "path not writeable" warnings when trying to update packages
  ## above the user's file permissions
  ## Use withCallingHandlers function to exit in case a package throws a warning,
  ## as suggested in https://stackoverflow.com/questions/26244530/how-do-i-make-install-packages-return-an-error-if-an-r-package-cannot-be-install
  withCallingHandlers(biocLite(not_installed_packages, suppressUpdates = T),
                      warning = function(w) stop(w))

  ## Message on successfull installs
  outmess <- c(outmess,"check_R_packages.R: Missing packages have been installed")
}

devtools::session_info()

for (i in outmess) {
  message(i)
}
