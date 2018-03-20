## Usage:
# Rscript gencode_gtf.R > gencode_gtf_log.txt 2>&1

library('rtracklayer')

## Get data from Gencode

## PRI regions
pri <- import('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.primary_assembly.annotation.gtf.gz')
length(pri)
pri_gene <- pri[pri$type == 'gene']
length(pri_gene)

## ALL regions
gall <- import('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz')
length(gall)
gall_gene <- gall[gall$type == 'gene']
length(gall_gene)


## Get an hg38 geneMap object
load('/dcl01/lieber/ajaffe/lab/libd_alzheimers/hg38_newidx/gene_info.Rdata')
nrow(geneMap)

## Compare against PRI
length(pri_gene) - nrow(geneMap)
sum(table(seqnames(pri_gene))[-(1:25)])
table(seqnames(pri_gene))[-(1:25)]

m <- match(geneMap$gencodeID, pri_gene$gene_id)
identical(m, seq_len(nrow(geneMap)))

## Check those in PRI but not in CHR
pri_gene[(nrow(geneMap) + 1):length(pri_gene)]
browseURL('http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000278704;r=GL000009.2:56140-58376;t=ENST00000618686')

rev_m <- match(pri_gene[(nrow(geneMap) + 1):length(pri_gene)]$gene_name, geneMap$Symbol)
table(is.na(rev_m))

## Check which ones are reverse matching
pri_gene[(nrow(geneMap) + 1):length(pri_gene)][!is.na(rev_m)]
geneMap[rev_m[!is.na(rev_m)], ]
browseURL('http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000274917;r=GL000220.1:112025-112177;t=ENST00000611446')
browseURL('http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000268674;r=KI270713.1:35407-35916;t=ENST00000601199')

## Compare against ALL
length(gall_gene) - nrow(geneMap)
sum(table(seqnames(gall_gene))[-(1:25)])
table(seqnames(gall_gene))[-(1:25)]

m2 <- match(geneMap$gencodeID, gall_gene$gene_id)
identical(m2, seq_len(nrow(geneMap)))

## Check those in ALL but not in CHR
gall_gene[(nrow(geneMap) + 1):length(gall_gene)]

rev_m2 <- match(gall_gene[(nrow(geneMap) + 1):length(gall_gene)]$gene_name, geneMap$Symbol)
table(is.na(rev_m2))
head(geneMap[rev_m2[!is.na(rev_m2)], ])

## Add number of times a gene in CHR is in PRI and ALL
geneCHR <- geneMap[, -which(colnames(geneMap) %in% c('meanExprs'))]
geneCHR$PRI_n <- 1
geneCHR$PRI_n[unique(rev_m[!is.na(rev_m)])] <- sapply(
    unique(rev_m[!is.na(rev_m)]), function(x) { 
        sum(x == rev_m[!is.na(rev_m)]) + 1 })

geneCHR$ALL_n <- 1
geneCHR$ALL_n[unique(rev_m2[!is.na(rev_m2)])] <- sapply(
    unique(rev_m2[!is.na(rev_m2)]), function(x) { 
        sum(x == rev_m2[!is.na(rev_m2)]) + 1 })
        
save(geneCHR, file = '/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE_PRI/geneCHR.Rdata')

table(geneCHR$PRI_n[geneCHR$PRI_n > 1])
sum(table(geneCHR$PRI_n[geneCHR$PRI_n > 1]))
table(geneCHR$ALL_n[geneCHR$ALL_n > 1])
sum(table(geneCHR$ALL_n[geneCHR$ALL_n > 1]))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
