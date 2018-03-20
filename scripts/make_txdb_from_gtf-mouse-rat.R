#######

library(GenomicFeatures)
library(GenomicRanges)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(rtracklayer)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")
library('getopt')
library('devtools')


## Specify parameters
spec <- matrix(c(
	'reference', 'r', 2, 'character', 'Indicate which reference genome should be used: [mm10, rn6]',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

if (opt$reference=="mm10") {

    ################################################################
    ####  mm10 - mouse  ############################################
    ################################################################

    ## chromosome info
    chrInfo = getChromInfoFromUCSC("mm10")
    chrInfo$chrom = as.character(chrInfo$chrom)
    chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:19,"X","Y","M")),]
    chrInfo$isCircular = rep(c(FALSE, TRUE), c(21,1))
    si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="mm10"))

    ##########################
    ##### GENCODE ############

    ## read in GTF as GRanges
    gencode_vM11 = import(con = "gencode.vM11.annotation.gtf", format = "gtf")
    seqlevels(gencode_vM11,force=TRUE) = paste0("chr", c(1:19,"X","Y","M"))
    seqinfo(gencode_vM11) = si

    # get map
    ensMap = mcols(gencode_vM11)
    ensMap = ensMap[!duplicated(ensMap$transcript_id),
	    c("gene_id","gene_name","transcript_id","gene_type")]

    ##  convert to txdb
    gencode_vM11_txdb = makeTxDbFromGRanges(gencode_vM11)

    # get introns
    introns = intronsByTranscript(gencode_vM11_txdb,use.names=TRUE)
    introns = unlist(introns)
    introns$TranscriptID = names(introns)
    introns$gencodeID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
    introns$ensemblID = ss(introns$gencodeID, "\\.")
    introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

    ### make unique junctions
    gencode_vM11_junctions = unique(introns)
    oo = findOverlaps(gencode_vM11_junctions, introns,type="equal")
    xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
    gencode_vM11_junctions$tx = CharacterList(xx)
    gencode_vM11_junctions$TranscriptID = NULL
    gencode_vM11_junctions$numTx = elementNROWS(gencode_vM11_junctions$tx)
    names(gencode_vM11_junctions) = NULL
    theJunctions = gencode_vM11_junctions
    save(theJunctions, file="junction_annotation_mm10_gencode_vM11.rda")


    ##########################
    ##### ENSEMBL ############

    ## read in GTF as GRanges
    ensembl_v86 = import(con = "Mus_musculus.GRCm38.86.gtf", format = "gtf")
    seqlevels(ensembl_v86,force=TRUE) = c(1:19,"X","Y","MT")
    seqlevels(ensembl_v86) = paste0("chr", c(1:19,"X","Y","M"))
    seqinfo(ensembl_v86) = si

    ##  convert to txdb
    ensembl_v86_txdb = makeTxDbFromGRanges(ensembl_v86)

    # get map
    ensMap = mcols(ensembl_v86)
    ensMap = ensMap[!duplicated(ensMap$transcript_id),
	    c("gene_id","gene_name","transcript_id","gene_biotype")]

    # get introns
    introns = intronsByTranscript(ensembl_v86_txdb,use.names=TRUE)
    introns = unlist(introns)
    introns$TranscriptID = names(introns)
    introns$ensemblID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
    introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

    ### make unique junctions
    ensembl_v86_junctions = unique(introns)
    oo = findOverlaps(ensembl_v86_junctions, introns,type="equal")
    xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
    ensembl_v86_junctions$tx = CharacterList(xx)
    ensembl_v86_junctions$TranscriptID = NULL
    ensembl_v86_junctions$numTx = elementNROWS(ensembl_v86_junctions$tx)
    names(ensembl_v86_junctions) = NULL
    theJunctions = ensembl_v86_junctions
    save(theJunctions, file="junction_annotation_mm10_ensembl_v86.rda")
}
if (opt$reference=="rn6") {

    ###################################################################
    ####  rn6 - rat  ##################################################
    ###################################################################

    ## chromosome info
    chrInfo = getChromInfoFromUCSC("rn6")
    chrInfo$chrom = as.character(chrInfo$chrom)
    chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:20,"X","Y","M")),]
    chrInfo$isCircular = rep(c(FALSE, TRUE), c(22,1))
    si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="rn6"))

    ##########################
    ##### ENSEMBL ############

    ## read in GTF as GRanges
    ensembl_v86 = import(con = "Rattus_norvegicus.Rnor_6.0.86.gtf", format = "gtf")
    seqlevels(ensembl_v86,force=TRUE) = c(1:20,"X","Y","MT")
    seqlevels(ensembl_v86) = paste0("chr", c(1:20,"X","Y","M"))
    seqinfo(ensembl_v86) = si

    ##  convert to txdb
    ensembl_v86_txdb = makeTxDbFromGRanges(ensembl_v86)

    # get map
    ensMap = mcols(ensembl_v86)
    ensMap = ensMap[!duplicated(ensMap$transcript_id),
	    c("gene_id","gene_name","transcript_id","gene_biotype")]

    # get introns
    introns = intronsByTranscript(ensembl_v86_txdb,use.names=TRUE)
    introns = unlist(introns)
    introns$TranscriptID = names(introns)
    introns$ensemblID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
    introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

    ### make unique junctions
    ensembl_v86_junctions = unique(introns)
    oo = findOverlaps(ensembl_v86_junctions, introns,type="equal")
    xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
    ensembl_v86_junctions$tx = CharacterList(xx)
    ensembl_v86_junctions$TranscriptID = NULL
    ensembl_v86_junctions$numTx = elementNROWS(ensembl_v86_junctions$tx)
    names(ensembl_v86_junctions) = NULL
    theJunctions = ensembl_v86_junctions
    save(theJunctions, file="junction_annotation_rn6_ensembl_v86.rda")
}

