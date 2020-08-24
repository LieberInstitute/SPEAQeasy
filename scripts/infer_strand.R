#  If R was installed locally, use packages associated with that installation
if (grepl("SPEAQeasy/Software/R-3.6.1/library", .libPaths()[1])) {
    repo_dir = dirname(dirname(dirname(.libPaths()[1])))
    library('checkpoint')
    checkpoint::checkpoint("2019-08-05",
                           project=paste0(repo_dir, '/scripts/'),
                           checkpointLocation = paste0(repo_dir, '/Software/R-3.6.1/'))
}

library('getopt')

#  These must sum to <= 0.5
thres_strand = 0.2   # How close to 1 the larger ratio (reverse or forward) must be to conclude strandedness
thres_unstrand = 0.1 # How close to 0.5 either ratio must be to be considered unstranded

spec <- matrix(c(
    'forward', 'f', 1, 'numeric', 'forward count of psuedo-aligned reads',
    'reverse', 'r', 1, 'numeric', 'forward count of psuedo-aligned reads',
    'paired', 'p', 1, 'character', '"paired" or "single"',
    'supposedStrand', 's', 1, 'character', 'user-provided strandness',
    'forceStrand', 'x', 1, 'logical', 'bypass errors from strand disagreement',
    'id', 'i', 1, 'character', 'sample ID'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

f_frac = opt$forward / (opt$forward + opt$reverse)
r_frac = opt$reverse / (opt$forward + opt$reverse)

#######################################################################
#  Determine the apparent strandness
#######################################################################

if (1 - f_frac < thres_strand) {
    strandness = "forward"
    print(paste0("Greater than ", round(100 * (1 - thres_strand), 2), "% of pseudoaligned reads matches forward-strandness."))
} else if (1 - r_frac < thres_strand) {
    strandness = "reverse"
    print(paste0("Greater than ", round(100 * (1 - thres_strand), 2), "% of pseudoaligned reads matches reverse-strandness."))
} else if (abs(f_frac - 0.5) < thres_unstrand) {
    strandness = "unstranded"
    print(paste0(round(100 * f_frac, 2), "% of pseudoaligned reads match forward-strandness;"))
    print(paste0(round(100 * r_frac, 2), "% of pseudoaligned reads match reverse-strandness."))
    print(paste0("The larger percentage is no more than ", round(100 * (0.5 + thres_unstrand), 2), "%, thus the reads appear unstranded."))
} else {
    default_message = paste0("based on psuedo-alignment testing, reads do not conclusively appear stranded nor unstranded.\n",
                             "          Observed ratios: ", round(f_frac, 3), " forward; ", round(r_frac, 3), " reverse. Thresholds used:\n",
                             "              Larger ratio > ", round(100 * (1 - thres_strand), 2), "% constitutes stranded reads;\n",
                             "              Ratios in [", round(100 * (0.5 - thres_unstrand ), 2), ", ", round(100 * (0.5 + thres_unstrand ), 2),
                             "]% constitute unstrandedness.\n")
    if (opt$forceStrand) {
        strandness = "unstranded"
        writeLines(paste0("Warning: ", default_message, " Pipeline execution will continue, as '--force_strand' was specified. ",
                          "This sample will be considered unstranded in the remaining analysis steps."))
    } else {
        writeLines(paste0(default_message, " This is a fatal error by default. To disable errors of this kind, ",
                          "you may add the command-line flag '--force_strand' and resume pipeline execution."))
        stop()
    }
}

#######################################################################
#  Throw an error or warning if apparent strandness doesn't agree with
#  user-provided strandness
#######################################################################

writeLines(strandness, con=paste0(opt$id, "_strandness_pattern.txt"))

if (strandness != opt$supposedStrand) {
    default_message = paste0("You have specified that reads should be ", opt$supposedStrand, 
                 "-stranded, but inference suggests ", strandness, "-strandness.")
    if (opt$forceStrand) {
        print(paste0("Warning: ", default_message, " Pipeline will resume, since you have provided '--force_strand'."))
    } else {
        stop(paste0(default_message, " Consider checking your samples for potential issues. Alternatively, you may add the command-line flag '--force_strand' and resume pipeline execution."))
    }
}
