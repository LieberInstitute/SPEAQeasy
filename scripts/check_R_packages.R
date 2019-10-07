## Start an outmess object to store messages
outmess <- "Log: checking R packages\n=============="

#  Use BiocManager to manage installation of Bioconductor packages: this
#  replaces the deprecated BiocInstaller and biocLite
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos='http://cran.us.r-project.org')

# Set the vector with package names required by the pipeline
# jaffelab was removed from this install block since it is installed in the previous block
packages <- c('devtools','Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db',
              'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
              'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db',
              'BSgenome.Rnorvegicus.UCSC.rn6', 'bumphunter',
              'getopt', 'BiocParallel','derfinder', 'rafalib', 'SummarizedExperiment',
              'plyr', 'rtracklayer', 'RColorBrewer','LieberInstitute/jaffelab')

## Try to load the required packages and send a message for each package that fails to load
load_res <- sapply(gsub('.*/', '', packages), requireNamespace, quietly = TRUE)
if(any(!load_res)) {
  outmess <- c(outmess,paste0("check_R_packages.R: missing the package: ",names(which(load_res == F))))
}

## Check if there are packages to be installed and procced accordingly
if( length(names(which(load_res == F))) == 0) {

  outmess <- c(outmess,"check_R_packages.R: Everything OK. No R packages need to be installed")

} else {

  outmess <- c(outmess,"check_R_packages.R: Attempting to install missing packages")

  ## Install missing packages: suppresing updates to avoid "path not writeable" warnings
  ## when trying to update packages above the user's file permissions
  present <- installed.packages(noCache = TRUE)
  BiocManager::install(packages[!gsub('.*/', '', packages) %in% rownames(present)],update=FALSE)
}

library('devtools') 
devtools::session_info()

for (i in outmess) {
  message(i)
}
