#  If R was installed locally, use packages associated with that installation
if (grepl("SPEAQeasy/Software/R-3.6.1/library", .libPaths()[1])) {
    repo_dir <- dirname(dirname(dirname(.libPaths()[1])))
    library("checkpoint")
    checkpoint::checkpoint("2019-08-05",
        project = paste0(repo_dir, "/scripts/"),
        checkpointLocation = paste0(repo_dir, "/Software/R-3.6.1/")
    )
}

library("Biostrings")
library("biomaRt")
library("jaffelab")
library("GenomicFeatures")
library("GenomicRanges")
library("org.Hs.eg.db")
library("rtracklayer")
library("devtools")
library("getopt")

#  The working directory should have the annotation .gtf file

spec <- matrix(c(
    "reference", "r", 1, "character", "hg38, hg19, mm10, or rn6",
    "suffix", "s", 1, "character", "suffix for filenames based on anno version"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

suffix <- paste0("_", opt$suffix)
temp <- strsplit(suffix, "_")[[1]]
opt$type <- temp[length(temp)] # "main", "primary", or "custom"

########################################################
#  Create chrom sizes file from assembly fasta
########################################################

fa <- readDNAStringSet(list.files(pattern = ".*\\.fa"))
fa_df <- data.frame(chr = names(fa), len = width(fa))
fa_df$chr <- as.character(fa_df$chr)
fa_df$chr <- ss(fa_df$chr, " ")

write.table(fa_df,
    file = paste0("chrom_sizes", suffix), row.names = FALSE,
    col.names = FALSE, quote = FALSE, sep = "\t"
)

########################################################
#  Create junction annotation object from gencode gtf
########################################################

if (opt$reference == "hg19" || opt$reference == "hg38") {
    chrom_names <- paste0("chr", c(1:22, "X", "Y", "M"))
} else if (opt$reference == "mm10") {
    chrom_names <- paste0("chr", c(1:19, "X", "Y", "M"))
} else { # rat/ rn6
    chrom_names <- paste0("chr", c(1:20, "X", "Y", "M"))
}

## chromosome info
chrInfo <- getChromInfoFromUCSC(opt$reference)
chrInfo$chrom <- as.character(chrInfo$chrom)

if (opt$type == "main") {
    chrInfo <- chrInfo[chrInfo$chrom %in% chrom_names, ]
    chrInfo$isCircular <- c(rep(FALSE, length(chrom_names) - 1), TRUE)
} else {
    isCircular <- rep(FALSE, length(unique(chrInfo$chrom)))
    isCircular[length(chrom_names)] <- TRUE
    chrInfo$isCircular <- isCircular
}
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome = opt$reference))

## read in GTF as GRanges
gencode_gtf <- import(con = list.files(pattern = ".*\\.gtf"), format = "gtf")
seqlevelsStyle(gencode_gtf) <- "UCSC"

seqlevels(gencode_gtf, pruning.mode = "coarse") <- seqlevels(si)
seqinfo(gencode_gtf) <- si

# get map
if (opt$reference == "rn6") {
    cols_to_select <- c("gene_id", "gene_name", "transcript_id")
} else {
    cols_to_select <- c("gene_id", "gene_name", "transcript_id", "gene_type")
}
ensMap <- mcols(gencode_gtf)
ensMap <- ensMap[!duplicated(ensMap$transcript_id), cols_to_select]

##  convert to txdb
gencode_txdb <- makeTxDbFromGRanges(gencode_gtf)

# get introns
introns <- intronsByTranscript(gencode_txdb, use.names = TRUE)
introns <- unlist(introns)
introns$TranscriptID <- names(introns)
if (opt$reference == "rn6") {
    introns$ensemblID <- ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
} else {
    introns$gencodeID <- ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
    introns$ensemblID <- ss(introns$gencodeID, "\\.")
}
introns$symbol <- ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

### make unique junctions
theJunctions <- unique(introns)
oo <- findOverlaps(theJunctions, introns, type = "equal")
xx <- split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
theJunctions$tx <- CharacterList(xx)
theJunctions$TranscriptID <- NULL
theJunctions$numTx <- elementNROWS(theJunctions$tx)
names(theJunctions) <- NULL
save(theJunctions, file = paste0("junction_annotation", suffix, ".rda"))

#####################################################################
#  Create feature-to-Tx object from gtf and junction annotation
#####################################################################

if (opt$reference != "rn6") {
    # get exons
    ensExons <- exonsBy(gencode_txdb)
    map <- select(gencode_txdb,
        keys = names(ensExons), keytype = "TXID",
        columns = c("TXID", "TXNAME", "GENEID")
    )
    ensExons <- unlist(ensExons)
    ensExons$TxID <- map$TXNAME[match(names(ensExons), map$TXID)]

    gtf_colnames <- c(
        "seqnames", "start", "end", "width", "strand", "gene_id",
        "gene_name"
    )
    geneMap <- data.frame(gencode_gtf)
    geneMap <- geneMap[geneMap$type == "gene", gtf_colnames]
    geneMap$ensemblID <- ss(geneMap$gene_id, "\\.")
    rownames(geneMap) <- geneMap$gene_id

    ## add gene
    geneToTx <- CharacterList(split(map$TXNAME, map$GENEID))
    geneMapGR <- makeGRangesFromDataFrame(geneMap, keep = TRUE)
    geneMapGR$GenTx <- CharacterList(vector("list", length(geneMapGR)))
    mmGene <- match(geneMapGR$gene_id, names(geneToTx))
    geneMapGR$GenTx[!is.na(mmGene)] <- geneToTx[mmGene[!is.na(mmGene)]]

    exonMap <- data.frame(gencode_gtf)
    exonMap <- exonMap[exonMap$type == "exon", gtf_colnames]
    exonMap$ensemblID <- ss(exonMap$gene_id, "\\.")
    eMap <- GRanges(exonMap$seqnames, IRanges(exonMap$start, exonMap$end))
    keepIndex <- which(!duplicated(eMap))
    exonMap <- exonMap[keepIndex, ]

    ## add exon
    exonMapGR <- makeGRangesFromDataFrame(exonMap, keep = TRUE)
    ooExon <- findOverlaps(exonMapGR, ensExons, type = "equal")
    exonToTx <- CharacterList(split(
        ensExons$TxID[subjectHits(ooExon)],
        names(exonMapGR[queryHits(ooExon)])
    ))
    exonMapGR$GenTx <- CharacterList(vector("list", length(exonMapGR)))
    mm <- match(names(exonMapGR), names(exonToTx))
    exonMapGR$GenTx[!is.na(mm)] <- exonToTx[mm[!is.na(mm)]]

    ## and junction:
    names(theJunctions) <- paste0(
        seqnames(theJunctions),
        ":", start(theJunctions), "-", end(theJunctions), "(",
        strand(theJunctions), ")"
    )
    theJunctions$nonStrandName <- paste0(
        seqnames(theJunctions),
        ":", start(theJunctions), "-", end(theJunctions), "(*)"
    )

    ############################
    #### features to tx ########
    gn <- geneMapGR$GenTx
    names(gn) <- names(geneMapGR)
    exn <- exonMapGR$GenTx
    names(exn) <- names(exonMapGR)
    jxn <- jxn2 <- theJunctions$tx
    names(jxn) <- names(theJunctions)
    names(jxn2) <- theJunctions$nonStrandName
    allTx <- c(gn, exn, jxn, jxn2)

    save(allTx, file = paste0("feature_to_Tx", suffix, ".rda"))
}
