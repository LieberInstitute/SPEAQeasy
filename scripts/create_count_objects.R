## Required libraries
library("derfinder")
library("BiocParallel")
library("Biostrings")
library("GenomicRanges")
library("GenomicFeatures")
library("jaffelab")
library("getopt")
library("rafalib")
library("devtools")
library("SummarizedExperiment")
library("DelayedArray")
library("matrixStats")
library("plyr")
library("rtracklayer")

## Specify parameters
spec <- matrix(c(
    "organism", "o", 2, "character", "'hg38', 'hg19', 'mm10', or 'rat'",
    "experiment", "e", 1, "character", "Experiment",
    "prefix", "p", 1, "character", "Prefix",
    "paired", "l", 1, "logical", "Whether the reads are paired-end or not",
    "ercc", "c", 1, "logical", "Whether the reads include ERCC or not",
    "cores", "t", 1, "integer", "Number of cores to use",
    "stranded", "s", 1, "character", "Strandedness of the data: Either 'FALSE', 'forward' or 'reverse'",
    "salmon", "n", 1, "logical", "Whether to use salmon quants rather than kallisto",
    "star", "r", 1, "logical", "Whether STAR was used for alignment",
    "output", "u", 1, "character", "Output directory for SPEAQeasy",
    "qsva_tx", "q", 2, "character", "Filename for QSVA TX list",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

#  Reference-specific libraries and sequence names
if (opt$organism == "hg38") {
    library("org.Hs.eg.db")

    chr_names <- paste0("chr", c(1:22, "X", "Y"))
    mito_chr <- "chrM"
} else if (opt$organism == "hg19") {
    library("org.Hs.eg.db")

    chr_names <- paste0("chr", c(1:22, "X", "Y"))
    mito_chr <- "chrM"
} else if (opt$organism == "mm10") {
    library("org.Mm.eg.db")

    chr_names <- paste0("chr", c(1:19, "X", "Y"))
    mito_chr <- "chrM"
} else { # 'rat'
    library("org.Rn.eg.db")

    chr_names <- c(1:20, "X", "Y")
    mito_chr <- "MT"
}


#  The default "prefix" is empty (i.e. "")
if (nchar(opt$prefix) > 0) {
    EXPNAME <- paste0(opt$experiment, "_", opt$prefix)
} else {
    EXPNAME <- opt$experiment
}

## read in pheno
manifest <- read.table("samples_complete.manifest",
    sep = " ",
    header = FALSE, stringsAsFactors = FALSE
)
metrics <- data.frame(
    "SAMPLE_ID" = manifest[, ncol(manifest) - 1],
    "strandness" = manifest[, ncol(manifest)],
    stringsAsFactors = FALSE
)
metrics$SAMPLE_ID <- as.character(metrics$SAMPLE_ID)

N <- length(metrics$SAMPLE_ID)


###############################################################################
#  Parse FastQC logs and summaries
###############################################################################

fastqcdata <- c(
    "SeqLength", "percentGC", "phred1", "phred2", "phred3", "phred4",
    "phredGT30", "phredGT35", "Adapter1", "Adapter2", "Adapter3"
)
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

#  Possible columns we expect in FastQC summaries
fastqc_stat_names <- c(
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
    "Kmer Content"
)

#  For choosing the post-trimming FastQC files for any samples that
#  were trimmed. If 'only_check_trim', then return a logical instead,
#  indicating if a sample was trimmed.
get_file <- function(id, suffix, read, only_check_trim = FALSE) {
    if (file.exists(paste0(id, read, "_trimmed", suffix))) {
        if (only_check_trim) {
            return(TRUE)
        } else {
            return(paste0(id, read, "_trimmed", suffix))
        }
    } else {
        if (only_check_trim) {
            return(FALSE)
        } else {
            return(paste0(id, read, "_untrimmed", suffix))
        }
    }
}

#  Parse a vector of FastQC summaries, containing each sample in the order
#  present in 'metrics'. Return a matrix of statistics
parse_summary <- function(summary_paths, metrics) {
    temp_rows <- lapply(summary_paths, function(x) {
        temp <- readLines(x)
        vals <- ss(temp, "\t")
        names(vals) <- ss(temp, "\t", 2)

        stopifnot(all(names(vals) %in% fastqc_stat_names))
        vals <- vals[fastqc_stat_names]
    })

    actual_colnames <- gsub(" ", "_", tolower(fastqc_stat_names))
    stat_mat <- matrix(unlist(temp_rows),
        ncol = length(fastqc_stat_names),
        byrow = TRUE,
        dimnames = list(metrics$SAMPLE_ID, actual_colnames)
    )

    return(stat_mat)
}

#  Parse a vector of FastQC "data reports", containing each sample in the order
#  present in 'metrics'. Return a matrix of statistics
parse_data <- function(data_path, metrics) {
    R <- lapply(data_path, function(x) {
        scan(x,
            what = "character", sep = "\n",
            quiet = TRUE, strip = TRUE
        )
    })
    names(R) <- metrics$SAMPLE_ID

    ## Split list into sublists of metric categories
    zz <- lapply(R, function(x) splitAt(x, which(x == ">>END_MODULE") + 1))

    # sequence length
    seqlen <- sapply(zz, function(x) {
        index <- grep(">>Basic Statistics", x)
        stopifnot(ss(x[[index]][9], "\t") == "Sequence length")
        return(ss(x[[index]][9], "\t", 2))
    })

    # percent GC
    gcp <- sapply(zz, function(x) {
        index <- grep(">>Basic Statistics", x)
        stopifnot(ss(x[[index]][10], "\t") == "%GC")
        return(ss(x[[index]][10], "\t", 2))
    })

    # median phred scores (at roughly 1/4, 1/2, 3/4, and end of seq length)
    # get positions
    index <- grep(">>Per base sequence quality", zz[[1]])
    len <- round((length(zz[[1]][[index]]) - 3) / 4)
    pos <- c(len + 3, 2 * len + 3, 3 * len + 3, length(zz[[1]][[index]]) - 1)
    nameSuf <- ss(zz[[1]][[index]][pos], "\t", 1)
    fastqcdata[3:6] <- paste0("phred", nameSuf)

    phred <- lapply(zz, function(x) {
        index <- grep(">>Per base sequence quality", x)
        stopifnot(ss(x[[index]][2], "\t", 3) == "Median")
        return(ss(x[[index]][pos], "\t", 3))
    })
    phred <- matrix(unlist(phred), ncol = 4, byrow = TRUE)

    # proportion of reads above phred 30 and 35
    sc <- lapply(zz, function(x) {
        return(x[[grep(">>Per sequence quality scores", x)]])
    })

    phred2 <- lapply(sc, function(x) x[3:(length(x) - 1)])
    phred2 <- lapply(phred2, function(x) {
        data.frame(score = ss(x, "\t", 1), count = ss(x, "\t", 2))
    })
    phred2 <- lapply(phred2, function(x) {
        data.frame(x, cumulRev = rev(cumsum(rev(as.numeric(levels(x$count))[x$count]))))
    })
    phred2 <- lapply(phred2, function(x) {
        data.frame(x, prop = x$cumulRev / x$cumulRev[1])
    })
    phred2 <- lapply(phred2, function(x) c(x[which(x$score %in% c(30, 35)), "prop"], 0, 0)[1:2])
    phred2 <- matrix(unlist(phred2), ncol = 2, byrow = TRUE)

    # Illumina adapter content (at roughly 1/2, 3/4, and end of seq length)
    # get positions
    ac <- lapply(zz, function(x) {
        return(x[[grep(">>Adapter Content", x)]])
    })

    len <- round((length(ac[[1]]) - 3) / 5)
    pos <- c(3 * len + 2, 4 * len + 2, length(ac[[1]]) - 1)
    nameSuf <- ss(ac[[1]][pos], "\t", 1)
    fastqcdata[9:11] <- paste0("Adapter", nameSuf)

    adap <- lapply(ac, function(x) x[pos])
    adap <- lapply(adap, function(x) ss(x, "\t", 2))
    adap <- matrix(as.numeric(unlist(adap)), ncol = 3, byrow = TRUE)

    combined <- data.frame(SeqLen = unlist(seqlen), GCprec = unlist(gcp), phred, phred2, adap)
    rownames(combined) <- NULL

    return(list(combined, fastqcdata))
}

#  Parse FastQC summaries and "data reports"; append to 'metrics'
if (opt$paired) {
    summary_paths_1 <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "_1")
    summary_paths_2 <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "_2")

    #  Get stat values for each mate
    stat_mat <- parse_summary(summary_paths_1, metrics)
    stat_mat_2 <- parse_summary(summary_paths_2, metrics)

    #  Combine into a single matrix and simplify some values
    stat_mat[] <- paste(stat_mat, stat_mat_2, sep = "/")
    stat_mat[stat_mat == "PASS/PASS"] <- "PASS"
    stat_mat[stat_mat == "WARN/WARN"] <- "WARN"
    stat_mat[stat_mat == "FAIL/FAIL"] <- "FAIL"
    stat_mat[stat_mat == "NA/NA"] <- "NA"

    metrics <- cbind(metrics, stat_mat)

    for (i in 1:2) {
        data_path <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_fastqc_data.txt", read = paste0("_", i))
        temp <- parse_data(data_path, metrics)
        combined <- temp[[1]]
        fastqcdata <- temp[[2]]
        names(combined) <- paste0(fastqcdata, "_R", i)
        metrics <- cbind(metrics, combined)
    }
} else {
    summary_paths <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "")
    metrics <- cbind(metrics, parse_summary(summary_paths, metrics))

    data_paths <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_fastqc_data.txt", read = "")
    temp <- parse_data(data_paths, metrics)
    combined <- temp[[1]]
    fastqcdata <- temp[[2]]
    rownames(combined) <- NULL
    names(combined) <- fastqcdata

    metrics <- cbind(metrics, combined)
}

#  Add a metric indicating if each sample was trimmed
if (opt$paired) {
    metrics$trimmed <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "_1", only_check_trim = TRUE)
} else {
    metrics$trimmed <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "", only_check_trim = TRUE)
}

sampIDs <- as.vector(metrics$SAMPLE_ID)

###############################################################################
#  Read in transcript pseudo-alignment stats, except for rat
###############################################################################

if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    if (opt$salmon) {
        #----------------------------------------------------------------------
        #  Salmon quantification
        #----------------------------------------------------------------------

        ## observed tpm and number of reads
        txTpm <- bplapply(sampIDs, function(x) {
            read.table(file.path(".", paste0(x, "_quant.sf")), header = TRUE)$TPM
        },
        BPPARAM = MulticoreParam(opt$cores)
        )
        txTpm <- do.call(cbind, txTpm)

        txNumReads <- bplapply(sampIDs, function(x) {
            read.table(file.path(".", paste0(x, "_quant.sf")), header = TRUE)$NumReads
        },
        BPPARAM = MulticoreParam(opt$cores)
        )
        txNumReads <- do.call(cbind, txNumReads)

        ## get names of transcripts
        txNames <- read.table(file.path(".", paste0(sampIDs[1], "_quant.sf")),
            header = TRUE
        )$Name
        txNames <- as.character(txNames)
    } else {
        #----------------------------------------------------------------------
        #  Kallisto quantification
        #----------------------------------------------------------------------

        txMatrices <- bplapply(sampIDs, function(x) {
            read.table(paste0(x, "_abundance.tsv"), header = TRUE)
        },
        BPPARAM = MulticoreParam(opt$cores)
        )

        txTpm <- do.call(cbind, lapply(txMatrices, function(x) x$tpm))
        txNumReads <- do.call(cbind, lapply(txMatrices, function(x) x$est_counts))
        txNames <- as.character(txMatrices[[1]]$target_id)
        rm(txMatrices)
    }

    colnames(txTpm) <- colnames(txNumReads) <- sampIDs

    txMap <- t(ss(txNames, "\\|", c(1, 7, 2, 6, 8)))
    txMap <- as.data.frame(txMap)
    rm(txNames)

    colnames(txMap) <- c("gencodeTx", "txLength", "gencodeID", "Symbol", "gene_type")
    rownames(txMap) <- rownames(txTpm) <- rownames(txNumReads) <- txMap$gencodeTx
}


###############################################################################
#  ERCC plots
###############################################################################

if (opt$ercc == TRUE) {

    ## observed kallisto tpm
    erccTPM <- sapply(sampIDs, function(x) {
        read.table(paste0(x, "_ercc_abundance.tsv"), header = TRUE)$tpm
    })
    rownames(erccTPM) <- read.table(paste0(sampIDs[1], "_ercc_abundance.tsv"),
        header = TRUE
    )$target_id
    # check finiteness / change NaNs to 0s
    erccTPM[which(is.na(erccTPM), arr.ind = T)] <- 0

    # expected concentration
    spikeIns <- read.delim("./ercc_actual_conc.txt",
        as.is = TRUE, row.names = 2
    )
    ## match row order
    spikeIns <- spikeIns[match(rownames(erccTPM), rownames(spikeIns)), ]

    pdf(file.path(".", "ercc_spikein_check_mix1.pdf"), h = 12, w = 18)
    mypar(4, 6)
    for (i in 1:ncol(erccTPM)) {
        plot(log2(10 * spikeIns[, "concentration.in.Mix.1..attomoles.ul."] + 1) ~ log2(erccTPM[, i] + 1),
            xlab = "Kallisto log2(TPM+1)", ylab = "Mix 1: log2(10*Concentration+1)",
            main = colnames(erccTPM)[i],
            xlim = c(min(log2(erccTPM + 1)), max(log2(erccTPM + 1)))
        )
        abline(0, 1, lty = 2)
    }
    dev.off()

    mix1conc <- matrix(rep(spikeIns[, "concentration.in.Mix.1..attomoles.ul."]),
        nc = ncol(erccTPM), nr = nrow(erccTPM), byrow = FALSE
    )
    logErr <- (log2(erccTPM + 1) - log2(10 * mix1conc + 1))
    metrics$ERCCsumLogErr <- colSums(logErr)
}

###############################################################################
#  Process alignment logs from HISAT2 or STAR, as applicable
###############################################################################

#  Return a desired value from a HISAT2 alignment summary.
#
#  Here 'log_lines' is a whitespace-stripped character vector containing each
#  line in the summary; 'phrase' is a string that should be present in the
#  lines of interest; 'get_perc' is a logical determining whether to extract
#  the percentage enclosed in parenthesis (rather than the first integer
#  present in the line); 'num_hits' is the number of times 'phrase' should be
#  present in a line.
parse_hisat_line <- function(log_lines, phrase, get_perc, num_hits = 1) {
    index <- grep(phrase, log_lines)
    if (length(index) != num_hits) {
        stop(
            paste0(
                "Phrase '", phrase, "' in HISAT2 log was expected ", num_hits,
                " time(s) but was observed ", length(index), " times."
            )
        )
    }

    if (get_perc) {
        #  Grab a percentage enclosed in parenthesis
        val <- ss(ss(log_lines[index], "\\(", 2), "%")
    } else {
        #  Grab a raw integer (the first number in every line)
        val <- gsub("%", "", ss(log_lines[index], " "))
    }

    val <- sum(as.numeric(val))

    return(val)
}

hisatStats <- function(logFile) {
    y <- scan(logFile,
        what = "character", sep = "\n",
        quiet = TRUE, strip = TRUE
    )
    if (opt$paired) {
        perc_paired <- parse_hisat_line(y, "%) were paired; of these:", TRUE)
        if (perc_paired == 100) {
            ## 100% of reads paired
            reads <- 2 * parse_hisat_line(y, "%) were paired; of these:", FALSE)
            unaligned <- parse_hisat_line(y, "%) aligned 0 times", FALSE)
        } else {
            ## Combo of paired and unpaired (from trimming w/ '--keep_unpaired')
            reads <- 2 * parse_hisat_line(y, "%) were paired; of these:", FALSE) +
                parse_hisat_line(y, "%) were unpaired; of these:", FALSE)
            unaligned <- parse_hisat_line(y, "%) aligned 0 times", FALSE, 2)
        }
    } else {
        reads <- parse_hisat_line(y, "reads; of these:", FALSE)
        unaligned <- parse_hisat_line(y, "%) aligned 0 times", FALSE)
    }

    o <- data.frame(
        numReads = reads,
        numMapped = reads - unaligned,
        numUnmapped = unaligned,
        overallMapRate = parse_hisat_line(y, "% overall alignment rate", FALSE) / 100
    )

    if (opt$paired) {
        o$concordMapRate <- (parse_hisat_line(y, "%) aligned concordantly exactly 1 time", TRUE) +
            parse_hisat_line(y, "%) aligned concordantly >1 times", TRUE)) / 100
    }

    return(o)
}

#  Return a 1-row data frame of STAR alignment metrics, given a single
#  sample ID and a character vector of metric names to extract. Note this
#  vector is "fixed" (the function is hardcoded to work for a particular
#  value of 'metric_names').
starStats <- function(id, metric_names) {
    #  Infer log path from sample ID, then read
    star_log <- readLines(paste0(id, "_STAR_alignment.log"))

    #  Infer whether trimming was performed from whether a post-trimming FastQC
    #  report exists for this sample ID
    if (opt$paired) {
        is_trimmed <- file.exists(paste0(id, "_1_trimmed_summary.txt"))
    } else {
        is_trimmed <- file.exists(paste0(id, "_trimmed_summary.txt"))
    }

    #  Grab the lines including each metric of interest
    key_lines <- star_log[sapply(metric_names, function(n) grep(n, star_log))]

    #  Extract the numeric value for each metric in those lines
    metric_values <- as.numeric(ss(key_lines, "\t", 2))

    o <- data.frame(
        "trimmed" = is_trimmed,
        "numReads" = metric_values[1],
        "numMapped" = sum(metric_values[5:6]),
        "numUnmapped" = sum(metric_values[2:4]),
        "overallMapRate" = sum(metric_values[5:6]) / metric_values[1]
    )

    return(o)
}

#  Extract alignment stats, and add to 'metrics'
if (opt$star) {
    metric_names <- c(
        "Number of input reads",
        "Number of reads unmapped: too many mismatches",
        "Number of reads unmapped: too short",
        "Number of reads unmapped: other",
        "Uniquely mapped reads number",
        "Number of reads mapped to multiple loci"
    )

    alignment_stats <- lapply(metrics$SAMPLE_ID, function(id) starStats(id, metric_names))
    alignment_stats <- do.call(rbind, alignment_stats)
} else { # HISAT2 is used for alignment
    logFiles <- file.path(".", paste0(metrics$SAMPLE_ID, "_align_summary.txt"))
    names(logFiles) <- metrics$SAMPLE_ID

    alignment_stats <- do.call(rbind, lapply(logFiles, hisatStats))
}

metrics <- cbind(metrics, alignment_stats)

###############################################################################
#  Process BAM files
###############################################################################

### confirm total mapping
bamFile <- paste0(metrics$SAMPLE_ID, "_sorted.bam")
metrics$bamFile <- file.path(opt$output, "alignment", "bam_sort", bamFile)

metrics$totalMapped <- unlist(bplapply(bamFile, getTotalMapped,
    chrs = chr_names,
    BPPARAM = MulticoreParam(opt$cores)
))
metrics$mitoMapped <- unlist(bplapply(bamFile, getTotalMapped,
    chrs = mito_chr,
    BPPARAM = MulticoreParam(opt$cores)
))
metrics$mitoRate <- metrics$mitoMapped / (metrics$mitoMapped + metrics$totalMapped)

###############################################################################
#  Read in GTF to get reference data regarding genes and exons
###############################################################################

gencodeGTF <- import(con = list.files(pattern = ".*\\.gtf"), format = "gtf")

gencodeGENES <- gencodeGTF[gencodeGTF$type == "gene", ]
names(gencodeGENES) <- gencodeGENES$gene_id

if (opt$organism == "mm10") {
    gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "exon_id")]
    names(gencodeEXONS) <- c("Chr", "Start", "End", "exon_gencodeID")
} else if (opt$organism == "rat") {
    gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "exon_id")]
    names(gencodeEXONS) <- c("Chr", "Start", "End", "exon_ensemblID")
} else { # human
    gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "gene_id", "exon_id")]
    names(gencodeEXONS) <- c("Chr", "Start", "End", "gene_id", "exon_gencodeID")
}

if (opt$organism %in% c("hg19", "hg38")) {
    ### exons in PAR regions
    par_y <- grep("PAR_Y", gencodeEXONS$gene_id)
    gencodeEXONS$exon_gencodeID[par_y] <- paste0(gencodeEXONS$exon_gencodeID[par_y], "_PAR_Y")
    gencodeEXONS[, -match("gene_id", colnames(gencodeEXONS))]
}


###############################################################################
#  Read in gene counts
###############################################################################

#  Get filenames for gene counts in sample order present in manifest
geneFn <- rep("", length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    geneFn[i] <- list.files(
        pattern = paste0(
            "^", metrics$SAMPLE_ID[i], "_.*_Genes\\.counts$"
        )
    )
}
names(geneFn) <- metrics$SAMPLE_ID

### read in annotation ##
geneMap <- read.delim(geneFn[1], skip = 1, as.is = TRUE)[
    , c("Geneid", "Length")
]

## organize gene map
indices <- match(geneMap$Geneid, gencodeGENES$gene_id)
if (any(is.na(indices))) {
    stop("Not all genes observed in FeatureCounts output are in GTF.")
}

#   Read in gene coordinates and strand from GTF
geneMap$Chr <- as.character(seqnames(gencodeGENES)[indices])
geneMap$Start <- start(gencodeGENES)[indices]
geneMap$End <- end(gencodeGENES)[indices]
geneMap$Strand <- as.character(strand(gencodeGENES)[indices])
rownames(geneMap) <- geneMap$Geneid

if (opt$organism == "rat") {
    geneMap$ensemblID <- geneMap$Geneid
} else {
    geneMap$gencodeID <- geneMap$Geneid
    geneMap$ensemblID <- ss(geneMap$Geneid, "\\.")
    geneMap$gene_type <- gencodeGENES$gene_type[indices]
}
geneMap$Geneid <- NULL

#  Get the 'Symbol' and 'EntrezID' columns
if (opt$organism %in% c("hg19", "hg38")) {
    geneMap$Symbol <- gencodeGENES$gene_name[indices]
    geneMap$EntrezID <- mapIds(org.Hs.eg.db, geneMap$Symbol, "ENTREZID", "SYMBOL")
} else if (opt$organism == "mm10") {
    temp <- gencodeGENES$gene_name[indices]
    geneMap$EntrezID <- mapIds(org.Mm.eg.db, temp, "ENTREZID", "SYMBOL")
    geneMap$Symbol <- mapIds(org.Mm.eg.db, temp, "MGI", "SYMBOL")
} else { # 'rat'
    geneMap$Symbol <- gencodeGENES$gene_name[indices]
    geneMap$EntrezID <- mapIds(org.Rn.eg.db, geneMap$Symbol, "ENTREZID", "SYMBOL")
}

## counts
geneCountList <- bplapply(geneFn,
    function(x) {
        cat(".")
        read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[, 1]
    },
    BPPARAM = MulticoreParam(opt$cores)
)
geneCounts <- do.call("cbind", geneCountList)
rownames(geneCounts) <- rownames(geneMap)
geneCounts <- geneCounts[, metrics$SAMPLE_ID] # put in order

# number of reads assigned
geneStatList <- lapply(paste0(geneFn, ".summary"),
    read.delim,
    row.names = 1
)
geneStats <- do.call("cbind", geneStatList)
colnames(geneStats) <- metrics$SAMPLE_ID
metrics$totalAssignedGene <- as.numeric(geneStats[1, ] / colSums(geneStats))

#  Add all the other stats from featureCounts at the gene level
geneStats_t <- t(geneStats)
colnames(geneStats_t) <- paste0("gene_", colnames(geneStats_t))
metrics <- cbind(metrics, geneStats_t)

#  rRNA rate: the 'gene_type' column does not exist in the rat GTF
if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    metrics$rRNA_rate <- colSums(geneCounts[which(geneMap$gene_type == "rRNA"), ]) / colSums(geneCounts)
}

# make RPKM
bg <- matrix(rep(as.numeric(geneStats["Assigned", ])),
    nc = nrow(metrics),
    nr = nrow(geneCounts), byrow = TRUE
)
widG <- matrix(rep(geneMap$Length),
    nr = nrow(geneCounts),
    nc = nrow(metrics), byrow = FALSE
)
geneRpkm <- geneCounts / (widG / 1000) / (bg / 1e6)

## save metrics
write.csv(
    metrics,
    file = paste0(
        "read_and_alignment_metrics_", EXPNAME, ".csv"
    )
)

###############################################################################
#  Read in exon counts
###############################################################################

#  Get filenames for exon counts in sample order present in manifest
exonFn <- rep("", length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    exonFn[i] <- list.files(
        pattern = paste0(
            "^", metrics$SAMPLE_ID[i], "_.*_Exons\\.counts$"
        )
    )
}
names(exonFn) <- metrics$SAMPLE_ID

### read in annotation ##
exonMap <- read.delim(exonFn[1], skip = 1, as.is = TRUE)[, 1:6]
rownames(exonMap) <- paste0("e", rownames(exonMap))

indices <- match(exonMap$Geneid, gencodeGENES$gene_id)
if (any(is.na(indices))) {
    stop("Not all genes of exons observed in FeatureCounts output are in GTF.")
}

if (opt$organism == "rat") {
    exonMap$ensemblID <- exonMap$Geneid
} else {
    exonMap$gencodeID <- exonMap$Geneid
    exonMap$ensemblID <- ss(exonMap$Geneid, "\\.")
    exonMap$gene_type <- gencodeGENES$gene_type[indices]
}
exonMap$Geneid <- NULL

#  Get the 'Symbol' and 'EntrezID' columns
if (opt$organism %in% c("hg19", "hg38")) {
    exonMap$Symbol <- gencodeGENES$gene_name[indices]
    exonMap$EntrezID <- mapIds(org.Hs.eg.db, exonMap$Symbol, "ENTREZID", "SYMBOL")
} else if (opt$organism == "mm10") {
    temp <- gencodeGENES$gene_name[indices]
    exonMap$EntrezID <- mapIds(org.Mm.eg.db, temp, "ENTREZID", "SYMBOL")
    exonMap$Symbol <- mapIds(org.Mm.eg.db, temp, "MGI", "SYMBOL")
} else { # 'rat'
    exonMap$Symbol <- gencodeGENES$gene_name[indices]
    exonMap$EntrezID <- mapIds(org.Rn.eg.db, exonMap$Symbol, "ENTREZID", "SYMBOL")
}

#  Add GENCODE exon ID
exonMap <- join(exonMap, gencodeEXONS, type = "left", match = "first")

## counts
exonCountList <- bplapply(exonFn,
    function(x) {
        cat(".")
        read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[, 1]
    },
    BPPARAM = MulticoreParam(opt$cores)
)
exonCounts <- do.call("cbind", exonCountList)
rownames(exonCounts) <- rownames(exonMap)
exonCounts <- exonCounts[, metrics$SAMPLE_ID] # put in order

## remove duplicated
eMap <- GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))

## drop runthrough exons with duplicated exons
i <- grepl("-", exonMap$Symbol)
j <- countOverlaps(eMap[i], eMap[!i], type = "equal") > 0
dropIndex <- which(i)[j]
if (length(dropIndex) > 0) {
    exonCounts <- exonCounts[-dropIndex, ]
    exonMap <- exonMap[-dropIndex, ]
    eMap <- eMap[-dropIndex, ]
}

## drop duplicated exons
keepIndex <- which(!duplicated(eMap))
exonCounts <- exonCounts[keepIndex, ]
exonMap <- exonMap[keepIndex, ]

#  Change rownames
exonMap$exon_libdID <- rownames(exonMap)

# number of reads assigned
exonStatList <- lapply(paste0(exonFn, ".summary"),
    read.delim,
    row.names = 1
)
exonStats <- do.call("cbind", exonStatList)
colnames(exonStats) <- metrics$SAMPLE_ID

## make RPKM
bgE <- matrix(rep(colSums(exonCounts)),
    nc = nrow(metrics),
    nr = nrow(exonCounts), byrow = TRUE
)
widE <- matrix(rep(exonMap$Length),
    nr = nrow(exonCounts),
    nc = nrow(metrics), byrow = FALSE
)
exonRpkm <- exonCounts / (widE / 1000) / (bgE / 1e6)

## add mean expression
geneMap$meanExprs <- rowMeans(geneRpkm)
exonMap$meanExprs <- rowMeans(exonRpkm)

###############################################################################
#  Add transcript maps
###############################################################################

if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    load(list.files(pattern = "feature_to_Tx.*\\.rda"))

    #  For gencode version 25 (the pipeline default), we have additional exon annotation
    if (file.exists("exonMaps_by_coord_hg38_gencode_v25.rda")) {
        load("exonMaps_by_coord_hg38_gencode_v25.rda")
    }

    ## gene annotation
    geneMap$Class <- "InGen"
    mmTx <- match(geneMap$gencodeID, names(allTx))
    tx <- CharacterList(vector("list", nrow(geneMap)))
    tx[!is.na(mmTx)] <- allTx[mmTx[!is.na(mmTx)]]
    geneMap$NumTx <- elementNROWS(tx)
    geneMap$gencodeTx <- sapply(tx, paste0, collapse = ";")

    ## exon annotation
    exonMap$Class <- "InGen"
    exonMap$coord <- paste0(exonMap$Chr, ":", exonMap$Start, "-", exonMap$End, "(", exonMap$Strand, ")")
    if (file.exists("exonMaps_by_coord_hg38_gencode_v25.rda")) {
        exonMap <- exonMap[, -which(colnames(exonMap) %in% c("exon_gencodeID", "exon_libdID"))]

        mmENSE <- match(exonMap$coord, names(coordToENSE))
        ENSE <- CharacterList(vector("list", nrow(exonMap)))
        ENSE[!is.na(mmENSE)] <- coordToENSE[mmENSE[!is.na(mmENSE)]]
        exonMap$NumENSE <- elementNROWS(ENSE)
        exonMap$exon_gencodeID <- sapply(ENSE, paste0, collapse = ";")

        mmLIBD <- match(exonMap$coord, names(coordToEid))
        libdID <- CharacterList(vector("list", nrow(exonMap)))
        libdID[!is.na(mmLIBD)] <- coordToEid[mmLIBD[!is.na(mmLIBD)]]
        exonMap$NumLIBD <- elementNROWS(libdID)
        exonMap$exon_libdID <- sapply(libdID, paste0, collapse = ";")

        mmTx <- match(exonMap$coord, names(coordToTX))
        tx <- CharacterList(vector("list", nrow(exonMap)))
        tx[!is.na(mmTx)] <- coordToTX[mmTx[!is.na(mmTx)]]
        exonMap$NumTx <- elementNROWS(tx)
        exonMap$gencodeTx <- sapply(tx, paste0, collapse = ";")
    } else {
        #  hg19 or hg38 without additional exon annotation for gencode release
        #  25, or mm10
        mmTx <- match(exonMap$gencodeID, names(allTx))
        tx <- CharacterList(vector("list", nrow(exonMap)))
        tx[!is.na(mmTx)] <- allTx[mmTx[!is.na(mmTx)]]
        exonMap$NumTx <- elementNROWS(tx)
        exonMap$gencodeTx <- sapply(tx, paste0, collapse = ";")
    }
}

###############################################################################
#  Create gene,exon RangedSummarizedExperiment objects
###############################################################################

#   SPEAQeasy settings were written to a CSV, which will be read in here and
#   used to populate the metadata of each RSE
meta_csv = read.csv('params.csv', header = FALSE)
rse_meta = meta_csv[,2]
names(rse_meta) = meta_csv[,1]
rse_meta = list('SPEAQeasy_settings' = as.list(rse_meta))

gr_genes <- GRanges(
    seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand
)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, -which(colnames(geneMap) %in%
    c("Chr", "Start", "End", "Strand"))])

rse_gene <- SummarizedExperiment(
    assays = list("counts" = geneCounts),
    rowRanges = gr_genes, colData = metrics, metadata = rse_meta
)
save(rse_gene, file = paste0("rse_gene_", EXPNAME, "_n", N, ".Rdata"))

gr_exons <- GRanges(
    seqnames = exonMap$Chr,
    IRanges(exonMap$Start, exonMap$End), strand = exonMap$Strand
)
names(gr_exons) <- rownames(exonMap)
mcols(gr_exons) <- DataFrame(exonMap[, -which(colnames(exonMap) %in%
    c("Chr", "Start", "End", "Strand"))])

rse_exon <- SummarizedExperiment(
    assays = list("counts" = exonCounts),
    rowRanges = gr_exons, colData = metrics, metadata = rse_meta
)
save(rse_exon, file = paste0("rse_exon_", EXPNAME, "_n", N, ".Rdata"))

###############################################################################
#  Junctions
###############################################################################

## import theJunctions annotation
load(list.files(pattern = "junction_annotation_.*\\.rda"))

## via primary alignments only
junctionFiles <- paste0(metrics$SAMPLE_ID, "_junctions_primaryOnly_regtools.count")

#  Use jaffelab::junctionCount. At this time, 'strandSpecific = TRUE' gives the
#  behavior we want even when data is unstranded
juncCounts <- junctionCount(junctionFiles, metrics$SAMPLE_ID,
    output = "Count", maxCores = opt$cores, strandSpecific = TRUE
)

## filter junction counts - drop jxns in <1% of samples
n <- max(1, floor(N / 100))
jCountsLogical <- DataFrame(sapply(juncCounts$countDF, function(x) x > 0))
jIndex <- which(rowSums(as.data.frame(jCountsLogical)) >= n)
juncCounts <- lapply(juncCounts, function(x) x[jIndex, ])

#------------------------------------------------------------------------------
#  anno/jMap
#------------------------------------------------------------------------------

anno <- juncCounts$anno

if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    ## add additional annotation
    anno$inGencode <- countOverlaps(anno, theJunctions, type = "equal") > 0
    anno$inGencodeStart <- countOverlaps(anno, theJunctions, type = "start") > 0
    anno$inGencodeEnd <- countOverlaps(anno, theJunctions, type = "end") > 0

    oo <- findOverlaps(anno, theJunctions, type = "equal")
    anno$gencodeGeneID <- NA
    anno$gencodeGeneID[queryHits(oo)] <- as.character(theJunctions$gencodeID[subjectHits(oo)])
    anno$ensemblID <- ss(anno$gencodeGeneID, "\\.")
    anno$Symbol <- NA
    anno$Symbol[queryHits(oo)] <- theJunctions$symbol[subjectHits(oo)]
    anno$gencodeStrand <- NA
    anno$gencodeStrand[queryHits(oo)] <- as.character(strand(theJunctions)[subjectHits(oo)])
    anno$gencodeTx <- CharacterList(vector("list", length(anno)))
    anno$gencodeTx[queryHits(oo)] <- theJunctions$tx[subjectHits(oo)]
    anno$numTx <- elementNROWS(anno$gencodeTx)

    ## junction code
    anno$code <- ifelse(anno$inGencode, "InGen",
        ifelse(anno$inGencodeStart & anno$inGencodeEnd, "ExonSkip",
            ifelse(anno$inGencodeStart | anno$inGencodeEnd, "AltStartEnd", "Novel")
        )
    )
} else { # 'rat'
    #  Rename to "UCSC"-style seq names
    seqlevelsStyle(anno) <- "UCSC"

    ## add additional annotation
    anno$inEnsembl <- countOverlaps(anno, theJunctions, type = "equal") > 0
    anno$inEnsemblStart <- countOverlaps(anno, theJunctions, type = "start") > 0
    anno$inEnsemblEnd <- countOverlaps(anno, theJunctions, type = "end") > 0

    oo <- findOverlaps(anno, theJunctions, type = "equal")
    anno$ensemblGeneID <- NA
    anno$ensemblGeneID[queryHits(oo)] <- as.character(theJunctions$ensemblID[subjectHits(oo)])
    anno$ensemblSymbol <- NA
    anno$ensemblSymbol[queryHits(oo)] <- theJunctions$symbol[subjectHits(oo)]
    anno$ensemblStrand <- NA
    anno$ensemblStrand[queryHits(oo)] <- as.character(strand(theJunctions)[subjectHits(oo)])
    anno$ensemblTx <- CharacterList(vector("list", length(anno)))
    anno$ensemblTx[queryHits(oo)] <- theJunctions$tx[subjectHits(oo)]
    anno$numTx <- elementNROWS(anno$ensemblTx)

    # clean up
    anno$ensemblSymbol <- geneMap$Symbol[match(anno$ensemblGeneID, rownames(geneMap))]

    ## junction code
    anno$code <- ifelse(anno$inEnsembl, "InEns",
        ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
            ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")
        )
    )
}

## b/w exons and junctions
exonGR <- GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
anno$startExon <- match(
    paste0(seqnames(anno), ":", start(anno) - 1),
    paste0(seqnames(exonGR), ":", end(exonGR))
)
anno$endExon <- match(
    paste0(seqnames(anno), ":", end(anno) + 1),
    paste0(seqnames(exonGR), ":", start(exonGR))
)

if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    g <- data.frame(
        leftGene = exonMap$gencodeID[anno$startExon],
        rightGene = exonMap$gencodeID[anno$endExon],
        leftGeneSym = exonMap$Symbol[anno$startExon],
        rightGeneSym = exonMap$Symbol[anno$endExon],
        stringsAsFactors = FALSE
    )
} else { # 'rat'
    g <- data.frame(
        leftGene = exonMap$ensemblID[anno$startExon],
        rightGene = exonMap$ensemblID[anno$endExon],
        leftGeneSym = exonMap$Symbol[anno$startExon],
        rightGeneSym = exonMap$Symbol[anno$endExon],
        stringsAsFactors = FALSE
    )
}

g$newGene <- NA
g$newGene[which(g$leftGene == g$rightGene)] <-
    g$leftGene[which(g$leftGene == g$rightGene)]
g$newGene[which(g$leftGene != g$rightGene)] <-
    paste0(g$leftGene, "-", g$rightGene)[which(g$leftGene != g$rightGene)]
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] <-
    g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))]
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] <-
    g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))]

anno$newGeneID <- g$newGene
anno$isFusion <- grepl("-", anno$newGeneID)

g$newGeneSym <- NA
g$newGeneSym[which(g$leftGene == g$rightGene)] <-
    g$leftGeneSym[which(g$leftGene == g$rightGene)]
g$newGeneSym[which(g$leftGene != g$rightGene)] <-
    paste0(g$leftGeneSym, "-", g$rightGeneSym)[which(g$leftGene != g$rightGene)]
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] <-
    g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))]
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] <-
    g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))]
g$newGeneSym[g$newGeneSym == ""] <- NA
g$newGeneSym[g$newGeneSym == "-"] <- NA

anno$newGeneSymbol <- g$newGeneSym

if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    anno$newGeneID[anno$code == "InGen"] <- anno$gencodeGeneID[anno$code == "InGen"]
    anno$newGeneSymbol[anno$code == "InGen"] <- anno$Symbol[anno$code == "InGen"]
}

## extract out jMap
jMap <- anno
colnames(mcols(jMap))[which(colnames(mcols(jMap)) == "code")] <- "Class"
rm(anno)

############ jCounts
jCounts <- as.matrix(as.data.frame(juncCounts$countDF))

#  Preserve exact colnames, which is actually not in general the case without
#  this line (columns/ sample IDs starting with an integer automatically get
#  an 'X' prepended, which we don't want)
colnames(jCounts) <- colnames(juncCounts$countDF)

jCounts <- jCounts[names(jMap), metrics$SAMPLE_ID]

############ jRpkm
bgJ <- matrix(rep(colSums(jCounts)),
    nc = nrow(metrics),
    nr = nrow(jCounts), byrow = TRUE
)
jRpkm <- jCounts / (bgJ / 10e6)

rownames(jCounts) <- rownames(jRpkm) <- names(jMap)
colnames(jRpkm) <- metrics$SAMPLE_ID
jMap$meanExprs <- rowMeans(jRpkm)

### save counts

tosaveCounts <- c("metrics", "geneCounts", "geneMap", "exonCounts", "exonMap", "jCounts", "jMap")
tosaveRpkm <- c("metrics", "geneRpkm", "geneMap", "exonRpkm", "exonMap", "jRpkm", "jMap")

#  Also save transcript data, except for rat
if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    tosaveCounts <- c(tosaveCounts, "txNumReads", "txMap")
    tosaveRpkm <- c(tosaveRpkm, "txNumReads", "txMap")
}

if (exists("erccTPM")) {
    tosaveCounts <- c("erccTPM", tosaveCounts)
    tosaveRpkm <- c("erccTPM", tosaveRpkm)
}

save(
    list = ls()[ls() %in% tosaveCounts], compress = TRUE,
    file = paste0("rawCounts_", EXPNAME, "_n", N, ".rda")
)
save(
    list = ls()[ls() %in% tosaveRpkm], compress = TRUE,
    file = paste0("rpkmCounts_", EXPNAME, "_n", N, ".rda")
)

rse_jx <- SummarizedExperiment(
    assays = list("counts" = jCounts), rowRanges = jMap, colData = metrics,
    metadata = rse_meta
)

save(rse_jx, file = paste0("rse_jx_", EXPNAME, "_n", N, ".Rdata"))

if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    ## transcript
    tx <- gencodeGTF[which(gencodeGTF$type == "transcript")]
    names(tx) <- tx$transcript_id

    #   Check that GTF contains observed transcripts
    if (!all((rownames(txTpm) %in% names(tx)))) {
        stop("Some transcripts do not appear to have corresponding annotation. If using custom annotation, please ensure the GTF has all transcripts present in the FASTA")
    }

    txMap <- tx[rownames(txTpm)]

    rse_tx <- SummarizedExperiment(
        assays = list("counts" = txNumReads, "tpm" = txTpm),
        colData = metrics, rowRanges = txMap, metadata = rse_meta
    )

    #  This file exists when the user specifies '--qsva'. Subset to
    #  user-specified transcripts in this case
    if (!is.null(opt$qsva_tx)) {
        select_tx <- readLines(opt$qsva_tx)
        if (!all(select_tx %in% rownames(rse_tx))) {
            stop("Selected transcripts passed via the '--qsva' argument are not all present in the final R object! Please check you are using appropriate transcript names (Ensembl ID), and check your annotation settings.")
        }

        rse_tx <- rse_tx[rownames(rse_tx) %in% select_tx, ]
    }

    save(rse_tx, file = paste0("rse_tx_", EXPNAME, "_n", N, ".Rdata"))
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
