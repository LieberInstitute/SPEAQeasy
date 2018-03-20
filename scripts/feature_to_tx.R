
###########################
###########################
###########################
############ hg19

### add transcript information
#library(GenomicFeatures)
#ensTxDb2 = loadDb("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/ensembl_v75_txdb.sqlite")
#seqlevels(ensTxDb2,force=TRUE) = c(1:22,"X","Y","MT")
#seqlevels(ensTxDb2) = paste0("chr", c(1:22,"X","Y","MT"))

### make txdb
library(GenomicFeatures)
library(org.Hs.eg.db)
library(rtracklayer)
library(biomaRt)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

## chromosome info
chrInfo = getChromInfoFromUCSC("hg19")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg19"))

##
ensembl_v75 = import(con = "../ensembl/Homo_sapiens.GRCh37.75.gtf", format = "gtf")
seqlevels(ensembl_v75,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(ensembl_v75) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(ensembl_v75) = si
ensTxDb = makeTxDbFromGRanges(ensembl_v75)
##
#gencode_v25 = import(con = "../GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf", format = "gtf")
#seqlevels(gencode_v25,force=TRUE) = c(1:22,"X","Y","M")
#seqlevels(gencode_v25) = paste0("chr", c(1:22,"X","Y","M"))
#seqinfo(gencode_v25) = si
#ensTxDb = makeTxDbFromGRanges(gencode_v25)

# get exons
ensExons = exonsBy(ensTxDb)
map = select(ensTxDb, keys=names(ensExons), keytype="TXID",
              columns=c("TXID", "TXNAME", "GENEID"))
ensExons = unlist(ensExons)
ensExons$TxID = map$TXNAME[match(names(ensExons),map$TXID)]


# load gene and exon maps
#load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
#load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/fullOct5/hg19/rpkmCounts_fulltest_oct5.hg19_n3.rda")

## create gene and exon maps from count files
geneFn = "./sample_counts/R16-099_Gencode.v25.hg19_Genes.counts"
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
geneMap$tempend = ""
geneMap$tempend[grep("_PAR_Y",geneMap$Geneid)] = "_PAR_Y"
geneMap$Geneid = paste0(ss(geneMap$Geneid, "\\."),geneMap$tempend)
rownames(geneMap) = geneMap$Geneid
geneMap$Geneid = NULL
geneMap$tempend = NULL

exonFn = "./sample_counts/R16-099_Gencode.v25.hg19_Exons.counts"
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$tempend = ""
exonMap$tempend[grep("_PAR_Y",exonMap$Geneid)] = "_PAR_Y"
exonMap$Geneid = paste0(ss(exonMap$Geneid, "\\."),exonMap$tempend)
exonMap$tempend = NULL

# biomaRt
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), values=rownames(geneMap), mart=ensembl)

geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]
exonMap$Symbol = sym$hgnc_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonMap = exonMap[keepIndex,]


## add gene
geneToTx = CharacterList(split(map$TXNAME, map$GENEID))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$EnsTx = CharacterList(vector("list", length(geneMapGR)))
mmGene = match(names(geneMapGR), names(geneToTx))
geneMapGR$EnsTx[!is.na(mmGene) ] = geneToTx[mmGene[!is.na(mmGene)]]

## add exon
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)	
ooExon = findOverlaps(exonMapGR, ensExons, type="equal")
exonToTx = CharacterList(split(ensExons$TxID[subjectHits(ooExon)],
	names(exonMapGR[queryHits(ooExon)])))
exonMapGR$EnsTx = CharacterList(vector("list", length(exonMapGR)))
mm = match(names(exonMapGR), names(exonToTx))
exonMapGR$EnsTx[!is.na(mm)] = exonToTx[mm[!is.na(mm)]]

## and junction:
#load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_hg19_v75_junction_annotation.rda")
load("junction_annotation_hg19_ensembl_v75.rda")

names(theJunctions) = paste0(seqnames(theJunctions),
	":", start(theJunctions), "-", end(theJunctions), "(",
	strand(theJunctions), ")")
theJunctions$nonStrandName = paste0(seqnames(theJunctions),
	":", start(theJunctions), "-", end(theJunctions), "(*)")

	
############################
#### features to tx ########
gn = geneMapGR$EnsTx
names(gn) = names(geneMapGR)
exn = exonMapGR$EnsTx
names(exn) = names(exonMapGR)
jxn = jxn2 = theJunctions$tx
names(jxn) = names(theJunctions)
names(jxn2) = theJunctions$nonStrandName
allTx = c(gn, exn, jxn,jxn2)

save(allTx, file="feature_to_Tx_ensembl_v75.rda")



###########################
###########################
###########################
############ hg38

### make txdb
library(GenomicFeatures)
library(org.Hs.eg.db)
library(rtracklayer)
library(biomaRt)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

## chromosome info
chrInfo = getChromInfoFromUCSC("hg38")
chrInfo$chrom = as.character(chrInfo$chrom)
chrInfo = chrInfo[chrInfo$chrom %in% paste0("chr", c(1:22,"X","Y", "M")),]
chrInfo$isCircular = rep(c(FALSE, TRUE), c(24,1))
si = with(chrInfo, Seqinfo(as.character(chrom), length, isCircular, genome="hg38"))

##
ensembl_v85 = import(con = "../ensembl/Homo_sapiens.GRCh38.85.gtf", format = "gtf")
seqlevels(ensembl_v85,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(ensembl_v85) = paste0("chr", c(1:22,"X","Y","M"))
seqinfo(ensembl_v85) = si
ensTxDb = makeTxDbFromGRanges(ensembl_v85)

# get exons
ensExons = exonsBy(ensTxDb)
map = select(ensTxDb, keys=names(ensExons), keytype="TXID",
              columns=c("TXID", "TXNAME", "GENEID"))
ensExons = unlist(ensExons)
ensExons$TxID = map$TXNAME[match(names(ensExons),map$TXID)]


## create gene and exon maps from count files
geneFn = "./sample_counts/R16-099_Gencode.v25.hg38_Genes.counts"
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
geneMap$tempend = ""
geneMap$tempend[grep("_PAR_Y",geneMap$Geneid)] = "_PAR_Y"
geneMap$Geneid = paste0(ss(geneMap$Geneid, "\\."),geneMap$tempend)
rownames(geneMap) = geneMap$Geneid
geneMap$Geneid = NULL
geneMap$tempend = NULL

exonFn = "./sample_counts/R16-099_Gencode.v25.hg38_Exons.counts"
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$tempend = ""
exonMap$tempend[grep("_PAR_Y",exonMap$Geneid)] = "_PAR_Y"
exonMap$Geneid = paste0(ss(exonMap$Geneid, "\\."),exonMap$tempend)
exonMap$tempend = NULL

# biomaRt
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), values=rownames(geneMap), mart=ensembl)

geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]
exonMap$Symbol = sym$hgnc_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonMap = exonMap[keepIndex,]


## add gene
geneToTx = CharacterList(split(map$TXNAME, map$GENEID))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$EnsTx = CharacterList(vector("list", length(geneMapGR)))
mmGene = match(names(geneMapGR), names(geneToTx))
geneMapGR$EnsTx[!is.na(mmGene) ] = geneToTx[mmGene[!is.na(mmGene)]]

## add exon
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)	
ooExon = findOverlaps(exonMapGR, ensExons, type="equal")
exonToTx = CharacterList(split(ensExons$TxID[subjectHits(ooExon)],
	names(exonMapGR[queryHits(ooExon)])))
exonMapGR$EnsTx = CharacterList(vector("list", length(exonMapGR)))
mm = match(names(exonMapGR), names(exonToTx))
exonMapGR$EnsTx[!is.na(mm)] = exonToTx[mm[!is.na(mm)]]

## and junction:
load("junction_annotation_hg38_ensembl_v85.rda")

names(theJunctions) = paste0(seqnames(theJunctions),
	":", start(theJunctions), "-", end(theJunctions), "(",
	strand(theJunctions), ")")
theJunctions$nonStrandName = paste0(seqnames(theJunctions),
	":", start(theJunctions), "-", end(theJunctions), "(*)")

	
############################
#### features to tx ########
gn = geneMapGR$EnsTx
names(gn) = names(geneMapGR)
exn = exonMapGR$EnsTx
names(exn) = names(exonMapGR)
jxn = jxn2 = theJunctions$tx
names(jxn) = names(theJunctions)
names(jxn2) = theJunctions$nonStrandName
allTx = c(gn, exn, jxn,jxn2)

save(allTx, file="feature_to_Tx_ensembl_v85.rda")
