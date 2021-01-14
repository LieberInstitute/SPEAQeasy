#  If R was installed locally, use packages associated with that installation
if (grepl("SPEAQeasy/Software/R-3.6.1/library", .libPaths()[1])) {
    repo_dir = dirname(dirname(dirname(.libPaths()[1])))
    library('checkpoint')
    checkpoint::checkpoint("2019-08-05",
                           project=paste0(repo_dir, '/scripts/'),
                           checkpointLocation = paste0(repo_dir, '/Software/R-3.6.1/'))
}

## Required libraries
library('derfinder')
library('BiocParallel')
library('jaffelab')
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'organism', 'o', 1, 'character', 'Either rn6, mm10 or human',
    'experiment', 'e', 1, 'character', 'Experiment',
    'cores', 'c', 1, 'integer', 'Number of cores to use',
    'strand', 's', 1, 'character', 'Strand of all bigwig files',
    'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## read in pheno
manifest <- read.table('samples_complete.manifest', sep = ' ', header = FALSE,
                       stringsAsFactors = FALSE)
ids = manifest[, ncol(manifest)-1]

#  Identify bigwig files
bigwigs <- sapply(ids, function(x) {
    list.files(pattern=paste0('^', x, '(|_(Reverse|Forward))\\.bw'))
})
stopifnot(length(bigwigs) == length(ids)) # one file per id is expected

## Chrs to use, mitocondrial chromosome has to be the last one for the code
## to work later on
if (opt$organism == "rn6") {
    CHR <- c(1:20,"X","Y","MT")
} else if (opt$organism == "mm10") {
    CHR <- paste0("chr",c(1:19,"X","Y","M"))
} else {
    CHR <- paste0("chr",c(1:22,"X","Y","M"))
}
    
#  Check bw file(s) to see which chromosomes have a nonzero number of mapped
#  reads. 'fullCoverage' will be passed this subset of chromosomes. If there
#  are multiple samples, the filtering process requires that all samples have a
#  >0 amount.
goodChrs = rep(TRUE, length(CHR))
for (f in bigwigs) {
    goodChrs = goodChrs &
               sapply(CHR, function(thisChr) getTotalMapped(f, thisChr) > 0)
}
CHR = CHR[goodChrs]

#  Perform "full coverage" from derfinder
fullCov <- fullCoverage(files = bigwigs, chrs = CHR,
                        mc.cores = opt$cores)

if (opt$strand == "unstranded") {
    filename = paste0('fullCoverage_', opt$experiment, '_n', length(ids),
                      '.rda')
} else {
    filename = paste0('fullCoverage_', opt$experiment, '_', opt$strand, '_n',
                      length(ids), '.rda')
}
save(fullCov, file = filename)


# Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
devtools::session_info()
