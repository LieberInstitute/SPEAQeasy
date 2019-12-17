library('getopt')

#  These must sum to <= 0.5
thres_strand = 0.2   # How close to 1 the larger ratio (reverse or forward) must be to conclude strandedness
thres_unstrand = 0.1 # How close to 0.5 either ratio must be to be considered unstranded

spec <- matrix(c(
    'forward', 'f', 1, 'numeric', 'forward count of psuedo-aligned reads',
    'reverse', 'r', 1, 'numeric', 'forward count of psuedo-aligned reads',
    'paired', 'p', 1, 'character', '"paired" or "single"',
    'supposedStrand', 's', 1, 'character', 'user-provided strandness'
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
    strandness = "unstranded"
    print("[Warning] Based on psuedo-alignment testing, reads do not conclusively appear stranded nor unstranded.")
    print(paste0("[Warning] Observed ratios: ", round(f_frac, 3), " forward; ", round(r_frac, 3), " reverse. Thresholds used:"))
    print(paste0("              Larger ratio > ", round(100 * (1 - thres_strand), 2), "% constitutes stranded reads;"))
    print(paste0("              Ratios in [", round(100 * (0.5 - thres_unstrand ), 2), ", ", round(100 * (0.5 + thres_unstrand ), 2),
                 "]% constitute unstrandedness."))
}

#######################################################################
#  Throw an error if apparent strandness is forward but user expects
#  reverse, or vice versa. Write strandness to file
#######################################################################

if (opt$paired == "paired") {
    id = strsplit(list.files(pattern=".+_1\\.f.*q.*"), "_1\\.f")[[1]][1]
} else {
    id = strsplit(list.files(pattern=".+\\.f.*q.*"), "\\.f")[[1]][1]
}

writeLines(strandness, con=paste0(id, "_strandness_pattern.txt"))

if (strandness != opt$supposedStrand && all(c(strandness, opt$supposedStrand) %in% c("reverse","forward"))) {
    stop(paste0("[Error] You have specified that reads should be ", opt$supposedStrand, 
                 "-stranded, but inference suggests ", strandness, 
                 "-strandness. Please re-check your samples or strand-specification before resuming."))
}
