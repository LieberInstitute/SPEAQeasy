library('Biostrings')
library('biomaRt')
library('jaffelab')
library('GenomicFeatures')
library('GenomicRanges')
library('org.Hs.eg.db')
library('rtracklayer')
library('devtools')
library('getopt')

#  The working directory should have the annotation .gtf file

spec <- matrix(c(
    'reference', 'r', 1, 'character', 'hg38, hg19, mm10, or rn6',
    'version', 'v', 1, 'character', 'release number for genome build',
    'type', 't', 1, 'character', '"primary" or "main" (referring to included chrs)'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

if (opt$reference == "hg19") {
    suffix = paste0('_', opt$reference, '_gencode_v', opt$version, 'lift37_', opt$type)
} else {
    suffix = paste0('_', opt$reference, '_gencode_v', opt$version, '_', opt$type)
}

########################################################
#  Create chrom sizes file from assembly fasta
########################################################

fa = readDNAStringSet(list.files(pattern=".*\\.fa"))
fa_df = data.frame(chr=names(fa), len=width(fa))
fa_df$chr = as.character(fa_df$chr)
fa_df$chr = ss(fa_df$chr, " ")

write.table(fa_df, file=paste0('chrom_sizes', suffix), row.names=FALSE, 
            col.names=FALSE, quote=FALSE, sep="\t")
            
########################################################
#  Create junction annotation object from gencode gtf
########################################################

if (opt$reference == "hg19" || opt$reference == "hg38") {
    chrom_names = paste0("chr", c(1:22,"X","Y", "M"))
} else if (opt$reference == "mm10") {
    chrom_names = paste0("chr", c(1:19,"X","Y","M"))
} else { # rat/ rn6
    chrom_names = paste0("chr", c(1:20, "X", "Y", "MT"))
}

## chromosome info
chrInfo = getChromInfoFromUCSC(opt$reference)
chrInfo$chrom = as.character(chrInfo$chrom)

if (opt$type == "main") {
    chrInfo = chrInfo[chrInfo$chrom %in% chrom_names,]
    chrInfo$isCircular = c(rep(FALSE, length(chrom_names)-1), TRUE) 
} else {
    isCircular = rep(FALSE, length(unique(chrInfo$chrom)))
    isCircular[length(chrom_names)] = TRUE
    chrInfo$isCircular = isCircular
}
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome=opt$reference))

## read in GTF as GRanges
gencode_gtf = import(con = list.files(pattern=".*\\.gtf"), format = "gtf")
seqlevelsStyle(gencode_gtf) = "UCSC"
seqinfo(gencode_gtf) = si

# get map
ensMap = mcols(gencode_gtf)
ensMap = ensMap[!duplicated(ensMap$transcript_id),
	    c("gene_id","gene_name","transcript_id","gene_type")]

##  convert to txdb
gencode_txdb = makeTxDbFromGRanges(gencode_gtf)

# get introns
introns = intronsByTranscript(gencode_txdb,use.names=TRUE)
introns = unlist(introns)
introns$TranscriptID = names(introns)
introns$gencodeID = ensMap$gene_id[match(introns$TranscriptID, ensMap$transcript_id)]
introns$ensemblID = ss(introns$gencodeID, "\\.")
introns$symbol = ensMap$gene_name[match(introns$TranscriptID, ensMap$transcript_id)]

### make unique junctions
gencode_junctions = unique(introns)
oo = findOverlaps(gencode_junctions, introns,type="equal")
xx = split(introns$TranscriptID[subjectHits(oo)], queryHits(oo))
gencode_junctions$tx = CharacterList(xx)
gencode_junctions$TranscriptID = NULL
gencode_junctions$numTx = elementNROWS(gencode_junctions$tx)
names(gencode_junctions) = NULL
#  Rename for compatibility with "create_count_objects-*.R" scripts
theJunctions = gencode_junctions
save(theJunctions, file=paste0("junction_annotation", suffix, ".rda"))

#####################################################################
#  Create feature-to-Tx object from gtf and junction annotation
#####################################################################

# get exons
ensExons = exonsBy(gencode_txdb)
map = select(gencode_txdb, keys=names(ensExons), keytype="TXID",
              columns=c("TXID", "TXNAME", "GENEID"))
ensExons = unlist(ensExons)
ensExons$TxID = map$TXNAME[match(names(ensExons),map$TXID)]


geneMap = data.frame(gencode_gtf)
geneMap = geneMap[geneMap$type == 'gene', c(1:5, 10, 13)]
geneMap$ensemblID = ss(geneMap$gene_id, "\\.")
rownames(geneMap) = geneMap$gene_id

## add gene
geneToTx = CharacterList(split(map$TXNAME, map$GENEID))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$GenTx = CharacterList(vector("list", length(geneMapGR)))
mmGene = match(geneMapGR$gene_id, names(geneToTx))
geneMapGR$GenTx[!is.na(mmGene) ] = geneToTx[mmGene[!is.na(mmGene)]]

exonMap = data.frame(gencode_gtf)
exonMap = exonMap[exonMap$type == 'exon', c(1:5, 10, 13)]
exonMap$ensemblID = ss(exonMap$gene_id, "\\.")
eMap = GRanges(exonMap$seqnames, IRanges(exonMap$start, exonMap$end))
keepIndex= which(!duplicated(eMap))
exonMap = exonMap[keepIndex,]

## add exon
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)	
ooExon = findOverlaps(exonMapGR, ensExons, type="equal")
exonToTx = CharacterList(split(ensExons$TxID[subjectHits(ooExon)],
	names(exonMapGR[queryHits(ooExon)])))
exonMapGR$GenTx = CharacterList(vector("list", length(exonMapGR)))
mm = match(names(exonMapGR), names(exonToTx))
exonMapGR$GenTx[!is.na(mm)] = exonToTx[mm[!is.na(mm)]]

## and junction:
names(gencode_junctions) = paste0(seqnames(gencode_junctions),
	":", start(gencode_junctions), "-", end(gencode_junctions), "(",
	strand(gencode_junctions), ")")
gencode_junctions$nonStrandName = paste0(seqnames(gencode_junctions),
	":", start(gencode_junctions), "-", end(gencode_junctions), "(*)")

############################
#### features to tx ########
gn = geneMapGR$GenTx
names(gn) = names(geneMapGR)
exn = exonMapGR$GenTx
names(exn) = names(exonMapGR)
jxn = jxn2 = gencode_junctions$tx
names(jxn) = names(gencode_junctions)
names(jxn2) = gencode_junctions$nonStrandName
allTx = c(gn, exn, jxn,jxn2)

save(allTx, file=paste0("feature_to_Tx", suffix, ".rda"))
