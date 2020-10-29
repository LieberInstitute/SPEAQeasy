## Required libraries
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'sampleids', 's', 2, 'character', 'path to samples.manifest file',
    'outdir', 'o', 1, 'character', 'Full path to directory where the merged fastq files will be saved to',
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
if(paired) system('touch .paired_end')

## Find the file extensions
files <- manifest[, 1]
extensions <- c('fastq.gz', 'fq.gz', 'fastq', 'fq')
patterns <- paste0(extensions, '$')
ext_found <- sapply(files, function(file) {
    extensions[names(unlist(sapply(patterns, grep, file))) == patterns]
})
if(any(is.na(ext_found))) {
    system('touch find_sample_error')
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}
extensions <- split(ext_found, manifest[, ncol(manifest)])

## Check that extensions are the same per group
if(any(sapply(extensions, function(x) { length(unique(x)) }) != 1)) {
    system('touch find_sample_error')
    stop("For each sample name, the extensions of the fastq files to be merged have to be the same")
}

## Is merging required?
merged <- length(unique(manifest[, ncol(manifest)])) == nrow(manifest)
if(!merged) {
    system('touch .requires_merging')
    message(paste0(Sys.time(), ' creating .samples_unmerged.manifest'))
    system(paste('mv', opt$sampleids, file.path(dirname(opt$sampleids),
        '.samples_unmerged.manifest')))

    message(paste(Sys.time(), 'creating the new samples.manifest file with the merged samples'))

    ## Split according to the sample names
    file_groups <- split(manifest, manifest[, ncol(manifest)])
    extensions <- sapply(extensions, '[', 1)
    
    if(paired) {
        new_manifest <- data.frame(
            file.path(opt$outdir, paste0(names(file_groups), '_1.', extensions)),
            rep(0, length(file_groups)),
            file.path(opt$outdir, paste0(names(file_groups), '_2.',
                extensions)),
            rep(0, length(file_groups)),
            names(file_groups), stringsAsFactors = FALSE
        )
    } else {
        new_manifest <- data.frame(
            file.path(opt$outdir, paste0(names(file_groups), '.', extensions)),
            rep(0, length(file_groups)),
            names(file_groups), stringsAsFactors = FALSE
        )
    }
    ## Make names short, in case you want to interactively check the new manifest
    colnames(new_manifest) <- paste0('V', seq_len(ncol(new_manifest)))

    write.table(new_manifest, file = opt$sampleids, row.names = FALSE,
        col.names = FALSE, quote = FALSE, sep = ' ')
}

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()