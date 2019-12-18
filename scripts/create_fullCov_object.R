## Required libraries
library('derfinder')
library('BiocParallel')
library('jaffelab')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'organism', 'o', 1, 'character', 'Either rn6, mm10 or human',
    'experiment', 'e', 1, 'character', 'Experiment',
    'paired', 'l', 1, 'logical', 'Whether the reads are paired-end or not',
    'fullcov', 'f', 1, 'logical', 'Whether to create the full coverage object or not',
    'cores', 'c', 1, 'integer', 'Number of cores to use',
    'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

if(opt$fullcov) {
    ## read in pheno
    manifest <- read.table('samples_complete.manifest',
        sep = ' ', header = FALSE, stringsAsFactors = FALSE)
    info <- data.frame('SAMPLE_ID' = manifest[, ncol(manifest)-1],
        stringsAsFactors = FALSE)
    N <- length(info$SAMPLE_ID)
    
    #  Infer the strand used for the experiment
    strandness = manifest[,ncol(manifest)]
    strand_counts = sapply(unique(strandness), function(x) length(which(strandness == x)))
    bad_indices = which((strand_counts > 0) & (strand_counts < max(strand_counts)))
    if (length(bad_indices) > 1) {
        print(paste0("[Warning] ", sum(strand_counts[bad_indices]), 
                     " samples had an inferred strand different from the majority strandness."))
    }
    experiment_strand = unique(strandness)[match(max(strand_counts), strand_counts)]

    ## add bigwig and bam files
    info$bamFile <- paste0(info$SAMPLE_ID, '_accepted_hits.sorted.bam')
    info$bwFile <- paste0(info$SAMPLE_ID, '.bw')

    ## Chrs to use, mitocondrial chromosome has to be the last one for the code
    ## to work later on
    if (opt$organism == "rn6") {
        CHR <- c(1:20,"X","Y","MT")
    } else if (opt$organism == "mm10") {
        CHR <- paste0("chr",c(1:19,"X","Y","M"))
    } else {
        CHR <- paste0("chr",c(1:22,"X","Y","M"))
    }
    
    #  Check bw file(s) to see which chromosomes have a nonzero number of mapped reads.
    #  'fullCoverage' will be passed this subset of chromosomes. If there are multiple
    #  samples, the filtering process requires that all samples have a >0 amount.
    goodChrs = rep(TRUE, length(CHR))
    for (f in info$bwFile) {
      goodChrs = goodChrs & sapply(CHR, function(thisChr) getTotalMapped(f, thisChr) > 0)
    }
    CHR = CHR[goodChrs]
      
    ###################################################################

    ## Uses BAM files if the bigwigs are strand specific
    if (experiment_strand == "unstranded") {
        fullCov <- fullCoverage(files = info$bwFile, chrs = CHR,
            mc.cores = opt$cores)
    } else {
        warning('Using the BAM files instead of the strand-specific BigWigs. You might want to run fullCoverage on the strand-specific BigWigs for your analysis purposes')
        fullCov <- fullCoverage(files = info$bamFile, chrs = CHR,
            mc.cores = opt$cores)
    }

    save(fullCov, file = paste0('fullCoverage_', opt$experiment, '_n', N, '.rda'))
}


# Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
gotDevtools <- requireNamespace('devtools', quietly = TRUE)
if(gotDevtools) {
    devtools::session_info()
} else {
    sessionInfo()
}