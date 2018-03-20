## Required libraries
library('getopt')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'outdir', 'o', 2, 'character', 'Path where the output of infer_experiment.py was saved. Defaults to HISAT2_out/infer_strandess',
    'pattern', 'p', 2, 'character', 'Name of the pattern file. Defaults to inferred_strandness_pattern.txt',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

if(is.null(opt$outdir)) opt <- c(opt, list('outdir' = c(getwd())))
if(is.null(opt$pattern)) opt <- c(opt, list('pattern' = 'inferred_strandness_pattern.txt'))


## In case this step crashes, assume that the data is not stranded
##write.table('none', file = opt$pattern, row.names = FALSE,
##    col.names = FALSE, quote = FALSE)

## Infer strandness
strandfiles <- dir(opt$outdir, pattern = 'txt', full.names = TRUE)
names(strandfiles) <- gsub('.txt', '', dir(opt$outdir, pattern = 'txt'))

inferred_strandness <- do.call(rbind, lapply(strandfiles, function(sf) {
    message(paste(Sys.time(), 'processing', sf))
    info <- readLines(sf)
    
    if(any(grepl('Unknown Data type', info))) {
        warning(paste('Unknown data type for', sf))
        return(data.frame(
            infer_library = 'Unknown',
            infer_frac_undetermined = NA,
            infer_pattern1 = NA,
            infer_frac_pattern1 = NA,
            infer_pattern2 = NA,
            infer_frac_pattern2 = NA,
            stringsAsFactors = FALSE
        ))
    }
    
    explained <- grep('explained', info)
    data.frame(
        infer_library = gsub('This is | Data', '', info[grep('This is', info)]),
        infer_frac_undetermined = as.numeric(gsub('.*: ', '',
            info[grep('determine', info)])),
        infer_pattern1 = gsub('".*', '', gsub('.*by "', '',
            info[explained[1]])),
        infer_frac_pattern1 = as.numeric(gsub('.*: ', '', info[explained[1]])),
        infer_pattern2 = gsub('".*', '', gsub('.*by "', '',
            info[explained[2]])),
        infer_frac_pattern2 = as.numeric(gsub('.*: ', '', info[explained[2]])),
        stringsAsFactors = FALSE
    ) 
}))

## Print some info
lapply(inferred_strandness[, -grep('frac', colnames(inferred_strandness))], table, useNA = 'ifany')
summary(inferred_strandness[, grep('frac', colnames(inferred_strandness))])

save(inferred_strandness , 
    file = file.path(opt$outdir, 'inferred_strandness.Rdata'))

## Explore visually the results
pdf(file.path(opt$outdir, 'inferred_strandness.pdf'))
plot(0.5,0.5, ylim = c(0,1), xlim = c(0,1), bty="n", xlab="", ylab="") ## Blank plot
polygon(x=c(0,1,1), y=c(1,0,1), col="gray90")  ## not possible
polygon(x=c(0,1,.5), y=c(0,0,.5), col="darkseagreen2")  ## pattern 1 (forward)
polygon(x=c(0,0,.5), y=c(0,1,.5), col="lightcoral")  ## pattern 2 (reverse)
polygon(x=c(0,0,.4,.6,.2), y=c(0,.2,.6,.4,0), col="moccasin")  ## unstranded
segments(0,1,1,0, col = 'red', lwd = 2)
par(new=TRUE) ## plot the points
with(inferred_strandness, plot(infer_frac_pattern1, infer_frac_pattern2,
    xlab = paste('Fraction of reads with pattern 1:',
    names(sort(table(inferred_strandness$infer_pattern1)))[1]),
    ylab = paste('Fraction of reads with pattern 2:',
    names(sort(table(inferred_strandness$infer_pattern2)))[1]),
    ylim = c(0,1), xlim = c(0,1), bty="n"
))
legend(.7,.95, c("Pattern 1", "Pattern 2", "None"),
    pch=15, col=c("darkseagreen2","lightcoral","moccasin"))

boxplot(inferred_strandness[, grep('frac', colnames(inferred_strandness))],
    ylim = c(0, 1), col = 'light blue', ylab = 'Fraction of reads')
dev.off()

## Determine which pattern to use
observed_diff = mean(inferred_strandness$infer_frac_pattern1 - inferred_strandness$infer_frac_pattern2)

pattern <- ifelse(observed_diff > 0.2, 
    names(sort(table(inferred_strandness$infer_pattern1)))[1],
    ifelse(observed_diff < -0.2,
        names(sort(table(inferred_strandness$infer_pattern2)))[1], 'none'))
        
message(paste(Sys.time(), 'will use the following pattern for bam2wig:',
    pattern))
write.table(pattern, file = opt$pattern, row.names = FALSE, col.names = FALSE,
    quote = FALSE)

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
session_info()
