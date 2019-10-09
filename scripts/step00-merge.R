## Required libraries
library('getopt')
library('BiocParallel')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'sampleids', 's', 1, 'character', 'Path to the samples.manifest file',
	'outdir', 'o', 1, 'character', 'Full path to directory where the merged fastq files will be saved to',
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

manifest <- read.table(opt$sampleids, sep = ' ', header = FALSE,
    stringsAsFactors = FALSE)

## Is the data paired end?
paired <- ncol(manifest) > 3

## Find the file extensions
files <- manifest[, 1]
extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
ext_found <- sapply(files, function(file) {
    extensions[names(unlist(sapply(patterns, grep, file))) == patterns]
})
if(any(is.na(ext_found))) {
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}
extensions <- split(ext_found, manifest[, ncol(manifest)])

## Check that extensions are the same per group
if(any(sapply(extensions, function(x) { length(unique(x)) }) != 1)) {
    stop("For each sample name, the extensions of the fastq files to be merged have to be the same")
}

## Create the output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## Split according to the sample names
file_groups <- split(manifest, manifest[, ncol(manifest)])
extensions <- sapply(extensions, '[', 1)

merge_files <- function(file_names, new_file) {
    message(paste(Sys.time(), 'creating', new_file))
    call <- paste('cat', paste(file_names, collapse = ' '), '>', new_file)
    print(call)
    system(call)
}

res <- bpmapply(function(common, new_name, extension) {
  if (paired) {
    merge_files(common[, 1],
        file.path(opt$outdir, paste0(new_name, '_1.', extension)))
    merge_files(common[, 3],
        file.path(opt$outdir, paste0(new_name, '_2.', extension)))
  } else {
    merge_files(common[, 1],
        file.path(opt$outdir, paste0(new_name, '.', extension)))
  }
}, file_groups, names(file_groups), extensions, 
    BPPARAM = MulticoreParam(opt$cores))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()