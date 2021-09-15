#  If R was installed locally, use packages associated with that installation
if (grepl("SPEAQeasy/Software/R-3.6.1/library", .libPaths()[1])) {
    repo_dir <- dirname(dirname(dirname(.libPaths()[1])))
    library("checkpoint")
    checkpoint::checkpoint("2019-08-05",
        project = paste0(repo_dir, "/scripts/"),
        checkpointLocation = paste0(repo_dir, "/Software/R-3.6.1/")
    )
}

library("getopt")

#  These must sum to <= 0.5
thres_strand <- 0.2 # How close to 1 the larger ratio (reverse or forward) must be to conclude strandedness
thres_unstrand <- 0.1 # How close to 0.5 either ratio must be to be considered unstranded

spec <- matrix(c(
    "forward", "f", 1, "numeric", "forward count of psuedo-aligned reads",
    "reverse", "r", 1, "numeric", "forward count of psuedo-aligned reads",
    "paired", "p", 1, "character", '"paired" or "single"',
    "supposedStrand", "s", 1, "character", "user-provided strandness",
    "strandMode", "m", 1, "character", "how to handle strand disagreement",
    "id", "i", 1, "character", "sample ID"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

f_frac <- opt$forward / (opt$forward + opt$reverse)
r_frac <- opt$reverse / (opt$forward + opt$reverse)

#######################################################################
#  Determine the apparent strandness
#######################################################################

if (1 - f_frac < thres_strand) {
    strandness <- "forward"
    print(paste0("Greater than ", round(100 * (1 - thres_strand), 2), "% of pseudoaligned reads matches forward-strandness."))
} else if (1 - r_frac < thres_strand) {
    strandness <- "reverse"
    print(paste0("Greater than ", round(100 * (1 - thres_strand), 2), "% of pseudoaligned reads matches reverse-strandness."))
} else if (abs(f_frac - 0.5) < thres_unstrand) {
    strandness <- "unstranded"
    print(paste0(round(100 * f_frac, 2), "% of pseudoaligned reads match forward-strandness;"))
    print(paste0(round(100 * r_frac, 2), "% of pseudoaligned reads match reverse-strandness."))
    print(paste0("The larger percentage is no more than ", round(100 * (0.5 + thres_unstrand), 2), "%, thus the reads appear unstranded."))
} else {
    default_message <- paste0(
        "based on psuedo-alignment testing, reads do not conclusively appear stranded nor unstranded.\n",
        "          Observed ratios: ", round(f_frac, 3), " forward; ", round(r_frac, 3), " reverse. Thresholds used:\n",
        "              Larger ratio > ", round(100 * (1 - thres_strand), 2), "% constitutes stranded reads;\n",
        "              Ratios in [", round(100 * (0.5 - thres_unstrand), 2), ", ", round(100 * (0.5 + thres_unstrand), 2),
        "]% constitute unstrandedness.\n"
    )
    
    #  This will later be changed to "unstranded"
    strandness <- "unknown"
    
    if (opt$strandMode == "accept") {
        writeLines(paste0(
            "Warning: ", default_message, "This sample will be considered ",
            "unstranded in the remaining analysis steps, since ",
            "'--strand_mode 'accept'' was specified. "
        ))
    } else if (opt$strandMode == "declare"){
        writeLines(paste0(
            "Warning: ", default_message, "This sample will be considered ",
            opt$supposedStrand, " in the remaining analysis steps, since ",
            "'--strand_mode 'declare'' was specified. "
        ))
    } else { # opt$strandMode = "strict"
        writeLines(default_message)
        stop("This is a fatal error by default. To disable errors of this kind, consider using '--strand_mode 'accept'' or '--strand_mode 'declare'' and resume pipeline execution.")
    }
}

#######################################################################
#  Handle potential disagreement between user-asserted and
#  SPEAQeasy-inferred strandness. Write appropriate strandness value
#  to a file.
#######################################################################

if (strandness != "unknown" && strandness != opt$supposedStrand) {
    default_message <- paste0(
        "You have specified that reads should be ", opt$supposedStrand,
        "-stranded, but inference suggests ", strandness, "-strandness."
    )
    if (opt$strandMode == "accept") {
        warning(paste0(default_message, " This sample will be considered '", strandness, "' since you have provided '--strand_mode 'accept''."))
        print(paste0("This sample will be considered '", strandness, "' since you have provided '--strand_mode 'accept''."))
    } else if (opt$strandMode == "declare"){
        warning(paste0(default_message, " This sample will be considered '", opt$supposedStrand, "' since you have provided '--strand_mode 'declare''."))
    } else { # opt$strandMode = "strict"
        stop(paste0(default_message, " Consider checking your samples for potential issues. Alternatively, you may consider using '--strand_mode 'accept'' or '--strand_mode 'declare'' and resume pipeline execution."))
    }
}

#  If strandness wasn't conclusively determined, we treat the sample as
#  unstranded
if (strandness == "unknown") {
    strandness = "unstranded"
}

if (opt$strandMode == "declare") {
    #  Use the declared strandness rather than the inferred one
    strandness = opt$supposedStrand
}

writeLines(strandness, con = paste0(opt$id, "_strandness_pattern.txt"))
