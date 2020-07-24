## Prepare BED files for http://rseqc.sourceforge.net/#infer-experiment-py
# module load R/3.3.x
# Rscript prep_bed.R -> prep_bed_log.txt 2>&1

library('rtracklayer')
library('devtools')
library('getopt')
## Specify parameters
spec <- matrix(c(
	'file', 'f', 2, 'character', 'Name of the GTF File to Read',
        'name', 'n', 2, 'character', 'Name of Reference Genome Selected at Runtime',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

gtfs <- file.path(c(opt$file))
names(gtfs) <- c(opt$name)

for(i in seq_along(gtfs)) {
    message(paste(Sys.time(), 'processing', gtfs[i]))
    gr <- import(gtfs[i])
    if('score' %in% colnames(mcols(gr))) {
        mcols(gr) <- mcols(gr)[, -which(colnames(mcols(gr)) == 'score')]
    }
    gr <- gr[gr$type == 'gene']
    message(paste(Sys.time(), 'exporting', paste0(names(gtfs)[i], '.bed')))
    export(gr, paste0(names(gtfs)[i], '.bed'), format = 'bed')
}

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
