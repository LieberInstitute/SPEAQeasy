## Usage:
# qrsh -l bluejay,mem_free=50G,h_vmem=55G,h_fsize=100G
# Rscript create_gs_gencode_v25_hg19.R > create_gs_gencode_v25_hg19_log.txt 2>&1

## Based on /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/gs/create_gs_gencode_v25.R


library('GenomicFeatures')
library('GenomicRanges')
library('derfinder')
library('rtracklayer')
library('devtools')

## Get the chromosome info for hg19
chrInfo <- getChromInfoFromUCSC('hg19')
chrInfo$chrom <- as.character(chrInfo$chrom)
chrInfo <- chrInfo[chrInfo$chrom %in% paste0('chr', c(1:22, 'X', 'Y', 'M')), ]
chrInfo$isCircular <- rep(c(FALSE, TRUE), c(24, 1))
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
    genome = 'hg19'))

## Load the gtf info
gencode_v25 <- import('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh37_hg19/gencode.v25lift37.annotation.gtf')

## Add the chromosome lengths
seqlevels(gencode_v25,force=TRUE) <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
seqinfo(gencode_v25) <- si

## Create the txdb object
gencode_v25_txdb_hg19 <- makeTxDbFromGRanges(gencode_v25)
save(gencode_v25_txdb_hg19, file = 'gencode_v25_txdb_hg19.Rdata')

## Create the genomic state object
gs_gencode_v25_hg19 <- makeGenomicState(gencode_v25_txdb_hg19,
        chrs = c(1:22, 'X', 'Y', 'M'))

## Save the genomic state object
save(gs_gencode_v25_hg19, file = 'gs_gencode_v25_hg19.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
