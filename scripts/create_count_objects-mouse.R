#  If R was installed locally, use packages associated with that installation
if (grepl("SPEAQeasy/Software/R-3.6.1/library", .libPaths()[1])) {
    repo_dir <- dirname(dirname(dirname(.libPaths()[1])))
    library("checkpoint")
    checkpoint::checkpoint("2019-08-05",
        project = paste0(repo_dir, "/scripts/"),
        checkpointLocation = paste0(repo_dir, "/Software/R-3.6.1/")
    )
}

## Required libraries
library("derfinder")
library("BiocParallel")
library("GenomicRanges")
library("GenomicFeatures")
library("org.Mm.eg.db")
library("BSgenome.Mmusculus.UCSC.mm10")
library("jaffelab")
library("getopt")
library("devtools")
library("SummarizedExperiment")
library("plyr")

## Specify parameters
spec <- matrix(c(
    "organism", "o", 2, "character", "mm10",
    "experiment", "e", 1, "character", "Experiment",
    "prefix", "p", 1, "character", "Prefix",
    "paired", "l", 1, "logical", "Whether the reads are paired-end or not",
    "ercc", "c", 1, "logical", "Whether the reads include ERCC or not",
    "cores", "t", 1, "integer", "Number of cores to use",
    "stranded", "s", 1, "character", "Strandedness of the data: Either 'FALSE', 'forward' or 'reverse'",
    "salmon", "n", 1, "logical", "Whether to use salmon quants rather than kallisto",
    "star", "r", 1, "logical", "Whether STAR was used for alignment",
    "help", "h", 0, "logical", "Display help"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

stopifnot(opt$stranded %in% c("FALSE", "forward", "reverse"))

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
    stringsAsFactors = FALSE
)
N <- length(metrics$SAMPLE_ID)

############################################################
###### FastQC results
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
#  were trimmed
get_file <- function(id, suffix, read) {
    if (file.exists(paste0(id, read, "_trimmed", suffix))) {
        return(paste0(id, read, "_trimmed", suffix))
    } else {
        return(paste0(id, read, "_untrimmed", suffix))
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

sampIDs <- as.vector(metrics$SAMPLE_ID)

if (opt$salmon) {
    ############################################################
    ###### salmon quantification

    ## observed tpm and number of reads
    txTpm <- bplapply(sampIDs, function(x) {
        read.table(file.path(paste0(x, "_quant.sf")), header = TRUE)$TPM
    },
    BPPARAM = MulticoreParam(opt$cores)
    )
    txTpm <- do.call(cbind, txTpm)

    txNumReads <- bplapply(sampIDs, function(x) {
        read.table(file.path(paste0(x, "_quant.sf")), header = TRUE)$NumReads
    },
    BPPARAM = MulticoreParam(opt$cores)
    )
    txNumReads <- do.call(cbind, txNumReads)

    # get names of transcripts
    txNames <- read.table(file.path(".", paste0(sampIDs[1], "_quant.sf")),
        header = TRUE
    )$Name
    txNames <- as.character(txNames)
} else {
    #######################################################################
    #  Kallisto quantification
    #######################################################################
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


############################################################
###### ercc plots
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

    pdf(file.path("ercc_spikein_check_mix1.pdf"), h = 12, w = 18)
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
############################################################


### add bam file
metrics$bamFile <- file.path(paste0(metrics$SAMPLE_ID, "_sorted.bam"))

#  Define functions for extracting metrics from HISAT2 logs
if (opt$paired == TRUE) {
    hisatStats <- function(logFile) {
        y <- scan(logFile,
            what = "character", sep = "\n",
            quiet = TRUE, strip = TRUE
        )

        if (as.numeric(ss(ss(y[2], "\\(", 2), "%")) == 100) {
            ## 100% of reads paired
            reads <- as.numeric(ss(y[1], " ")) * 2
            unaligned <- as.numeric(ss(y[12], " "))
            o <- data.frame(
                trimmed = FALSE,
                numReads = reads,
                numMapped = reads - unaligned,
                numUnmapped = unaligned,
                overallMapRate = as.numeric(ss(y[15], "\\%")) / 100,
                concordMapRate = (as.numeric(ss(ss(y[4], "\\(", 2), "%")) + as.numeric(ss(ss(y[5], "\\(", 2), "%"))) / 100,
                stringsAsFactors = FALSE
            )
        } else {
            ## Combo of paired and unpaired (from trimming)
            reads <- as.numeric(ss(y[2], " ")) * 2 + as.numeric(ss(y[15], " "))
            unaligned <- as.numeric(ss(y[12], " ")) + as.numeric(ss(y[16], " "))
            o <- data.frame(
                trimmed = TRUE,
                numReads = reads,
                numMapped = reads - unaligned,
                numUnmapped = unaligned,
                overallMapRate = as.numeric(ss(y[19], "\\%")) / 100,
                concordMapRate = (as.numeric(ss(ss(y[4], "\\(", 2), "%")) + as.numeric(ss(ss(y[5], "\\(", 2), "%"))) / 100,
                stringsAsFactors = FALSE
            )
        }
    }
} else {
    ## all reads unpaired
    hisatStats <- function(logFile) {
        y <- scan(logFile,
            what = "character", sep = "\n",
            quiet = TRUE, strip = TRUE
        )
        o <- data.frame(
            numReads = as.numeric(ss(y[1], " ")),
            numMapped = as.numeric(ss(y[1], " ")) - as.numeric(ss(y[3], " ")),
            numUnmapped = as.numeric(ss(y[3], " ")),
            overallMapRate = as.numeric(ss(y[6], "\\%")) / 100
        )
    }
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

### confirm total mapping
metrics$totalMapped <- unlist(bplapply(metrics$bamFile, getTotalMapped,
    chrs = paste0("chr", c(1:19, "X", "Y")),
    BPPARAM = MulticoreParam(opt$cores)
))
metrics$mitoMapped <- unlist(bplapply(metrics$bamFile, getTotalMapped,
    chrs = "chrM",
    BPPARAM = MulticoreParam(opt$cores)
))
metrics$mitoRate <- metrics$mitoMapped / (metrics$mitoMapped + metrics$totalMapped)

###################################################################

gencodeGTF <- import(con = list.files(pattern = ".*\\.gtf"), format = "gtf")
gencodeGENES <- mcols(gencodeGTF)[which(gencodeGTF$type == "gene"), c("gene_id", "type", "gene_type", "gene_name")]
rownames(gencodeGENES) <- gencodeGENES$gene_id

gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "exon_id")]
names(gencodeEXONS) <- c("Chr", "Start", "End", "exon_gencodeID")


###############
### gene counts

#  Get filenames for gene counts in sample order present in manifest
geneFn <- rep("", length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    geneFn[i] <- list.files(pattern = paste0(
        metrics$SAMPLE_ID[i], "_",
        opt$organism, "_.*_Genes\\.counts$"
    ))
}
names(geneFn) <- metrics$SAMPLE_ID

### read in annotation ##
geneMap <- read.delim(geneFn[1], skip = 1, as.is = TRUE)[, 1:6]

## organize gene map
geneMap$Chr <- ss(geneMap$Chr, ";")
geneMap$Start <- as.numeric(ss(geneMap$Start, ";"))
tmp <- strsplit(geneMap$End, ";")
geneMap$End <- as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand <- ss(geneMap$Strand, ";")
rownames(geneMap) <- geneMap$Geneid
geneMap$gencodeID <- geneMap$Geneid
geneMap$ensemblID <- ss(geneMap$Geneid, "\\.")
geneMap$Geneid <- NULL
geneMap$gene_type <- gencodeGENES[geneMap$gencodeID, "gene_type"]

temp <- gencodeGENES[geneMap$gencodeID, "gene_name"]
geneMap$EntrezID <- mapIds(org.Mm.eg.db, temp, "ENTREZID", "SYMBOL")
geneMap$Symbol <- mapIds(org.Mm.eg.db, temp, "MGI", "SYMBOL")


## counts
geneCountList <- mclapply(geneFn, function(x) {
    cat(".")
    read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[, 1]
}, mc.cores = opt$cores)
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
# rna Rate
metrics$rRNA_rate <- colSums(geneCounts[which(geneMap$gene_type == "rRNA"), ]) / colSums(geneCounts)


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
write.csv(metrics, file = paste0("read_and_alignment_metrics_", EXPNAME, ".csv"))


###############
### exon counts

#  Get filenames for exon counts in sample order present in manifest
exonFn <- rep("", length(metrics$SAMPLE_ID))
for (i in 1:length(metrics$SAMPLE_ID)) {
    exonFn[i] <- list.files(pattern = paste0(
        metrics$SAMPLE_ID[i], "_",
        opt$organism, "_.*_Exons\\.counts$"
    ))
}
names(exonFn) <- metrics$SAMPLE_ID

### read in annotation ##
exonMap <- read.delim(exonFn[1], skip = 1, as.is = TRUE)[, 1:6]
exonMap$gencodeID <- exonMap$Geneid
exonMap$ensemblID <- ss(exonMap$Geneid, "\\.")
rownames(exonMap) <- paste0("e", rownames(exonMap))
exonMap$Geneid <- NULL
exonMap$gene_type <- gencodeGENES[exonMap$gencodeID, "gene_type"]

temp <- gencodeGENES[exonMap$gencodeID, "gene_name"]
exonMap$EntrezID <- mapIds(org.Mm.eg.db, temp, "ENTREZID", "SYMBOL")
exonMap$Symbol <- mapIds(org.Mm.eg.db, temp, "MGI", "SYMBOL")

## add gencode exon id
exonMap <- join(exonMap, gencodeEXONS, type = "left", match = "first")

## counts
exonCountList <- mclapply(exonFn, function(x) {
    cat(".")
    read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[, 1]
}, mc.cores = opt$cores)
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

## change rownames
exonMap$exon_libdID <- rownames(exonMap)
# rownames(exonMap) = rownames(exonCounts) = exonMap$exon_gencodeID

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

## Create gene,exon RangedSummarizedExperiment objects

gr_genes <- GRanges(
    seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand
)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, -which(colnames(geneMap) %in%
    c("Chr", "Start", "End", "Strand"))])

rse_gene <- SummarizedExperiment(
    assays = list("counts" = geneCounts),
    rowRanges = gr_genes, colData = metrics
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
    rowRanges = gr_exons, colData = metrics
)
save(rse_exon, file = paste0("rse_exon_", EXPNAME, "_n", N, ".Rdata"))


###################
##### junctions

## via primary alignments only
junctionFiles <- file.path(paste0(metrics$SAMPLE_ID, "_junctions_primaryOnly_regtools.count"))
stopifnot(all(file.exists(junctionFiles))) #  TRUE

## annotate junctions
load(list.files(pattern = "junction_annotation_.*\\.rda"))

#  Handle strand in a consistent way regardless of differences between samples-
#  if any samples are determined to be unstranded, process all samples as if unstranded.
#  Otherwise, process samples as all stranded.
if (any(manifest[, ncol(manifest)] == "unstranded")) {
    juncCounts <- junctionCount(junctionFiles, metrics$SAMPLE_ID,
        output = "Count", maxCores = opt$cores, strandSpecific = FALSE
    )
} else {
    juncCounts <- junctionCount(junctionFiles, metrics$SAMPLE_ID,
        output = "Count", maxCores = opt$cores, strandSpecific = TRUE
    )
}
## filter junction counts - drop jxns in <1% of samples
n <- max(1, floor(N / 100))
jCountsLogical <- DataFrame(sapply(juncCounts$countDF, function(x) x > 0))
jIndex <- which(rowSums(as.data.frame(jCountsLogical)) >= n)
juncCounts <- lapply(juncCounts, function(x) x[jIndex, ])


############ anno/jMap
anno <- juncCounts$anno

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

g <- data.frame(
    leftGene = exonMap$gencodeID[anno$startExon],
    rightGene = exonMap$gencodeID[anno$endExon],
    leftGeneSym = exonMap$Symbol[anno$startExon],
    rightGeneSym = exonMap$Symbol[anno$endExon],
    stringsAsFactors = FALSE
)

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
anno$newGeneID[anno$code == "InGen"] <- anno$gencodeGeneID[anno$code == "InGen"]

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
anno$newGeneSymbol[anno$code == "InGen"] <- anno$Symbol[anno$code == "InGen"]

## extract out jMap
jMap <- anno
colnames(mcols(jMap))[which(colnames(mcols(jMap)) == "code")] <- "Class"
rm(anno)

############ jCounts
jCounts <- as.matrix(as.data.frame(juncCounts$countDF))
jCounts <- jCounts[names(jMap), gsub("-", ".", metrics$SAMPLE_ID)] # ensure lines up
colnames(jCounts) <- metrics$SAMPLE_ID # change from '.' to hyphens if needed

############ jRpkm
bgJ <- matrix(rep(colSums(jCounts)),
    nc = nrow(metrics),
    nr = nrow(jCounts), byrow = TRUE
)
jRpkm <- jCounts / (bgJ / 10e6)

rownames(jCounts) <- rownames(jRpkm) <- names(jMap)
colnames(jRpkm) <- metrics$SAMPLE_ID
jMap$meanExprs <- rowMeans(jRpkm)


# ## sequence of acceptor/donor sites
# left = right = jMap
# end(left) = start(left) +1
# start(right) = end(right) -1
# jMap$leftSeq  = getSeq(Hsapiens, left)
# jMap$rightSeq = getSeq(Hsapiens, right)



### save counts

tosaveCounts <- c(
    "metrics", "geneCounts", "geneMap", "exonCounts", "exonMap", "jCounts", "jMap",
    "txNumReads", "txMap"
)
tosaveRpkm <- c(
    "metrics", "geneRpkm", "geneMap", "exonRpkm", "exonMap", "jRpkm", "jMap",
    "txTpm", "txNumReads", "txMap"
)

if (exists("erccTPM")) {
    tosaveCounts <- c("erccTPM", tosaveCounts)
    tosaveRpkm <- c("erccTPM", tosaveRpkm)
}

save(
    list = ls()[ls() %in% tosaveCounts], compress = TRUE,
    file = file.path(paste0("rawCounts_", EXPNAME, "_n", N, ".rda"))
)
save(
    list = ls()[ls() %in% tosaveRpkm], compress = TRUE,
    file = file.path(paste0("rpkmCounts_", EXPNAME, "_n", N, ".rda"))
)


## Create RangedSummarizedExperiment objects
rse_jx <- SummarizedExperiment(
    assays = list("counts" = jCounts),
    rowRanges = jMap, colData = metrics
)
save(rse_jx, file = paste0("rse_jx_", EXPNAME, "_n", N, ".Rdata"))

## transcript
tx <- gencodeGTF[which(gencodeGTF$type == "transcript")]
names(tx) <- tx$transcript_id
txMap <- tx[rownames(txTpm)]
rse_tx <- SummarizedExperiment(
    assays = list("tpm" = txTpm),
    colData = metrics, rowRanges = txMap
)
save(rse_tx, file = paste0("rse_tx_", EXPNAME, "_n", N, ".Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
