## Usage:
# qrsh -l bluejay,mem_free=50G,h_vmem=55G,h_fsize=100G
# Rscript create_gs_gencode_v25.R > create_gs_gencode_v25_log.txt 2>&1

## Based on /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb/make_txdb_from_gtf.R

library('GenomicFeatures')
library('GenomicRanges')
library('derfinder')
library('rtracklayer')
library('devtools')

## Get the chromosome info for hg38
chrInfo <- getChromInfoFromUCSC('hg38')
chrInfo$chrom <- as.character(chrInfo$chrom)
chrInfo <- chrInfo[chrInfo$chrom %in% paste0('chr', c(1:22, 'X', 'Y', 'M')), ]
chrInfo$isCircular <- rep(c(FALSE, TRUE), c(24, 1))
si <- with(chrInfo, Seqinfo(as.character(chrom), length, isCircular,
    genome = 'hg38'))

## Load the gtf info
gencode_v25 <- import('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf')

## Add the chromosome lengths
seqlevels(gencode_v25,force=TRUE) <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
seqinfo(gencode_v25) <- si

## Create the txdb object
gencode_v25_txdb_hg38 <- makeTxDbFromGRanges(gencode_v25)
save(gencode_v25_txdb_hg38, file = 'gencode_v25_txdb_hg38.Rdata')

## Create the genomic state object
gs_gencode_v25_hg38 <- makeGenomicState(gencode_v25_txdb_hg38,
        chrs = c(1:22, 'X', 'Y', 'M'))

## Save the genomic state object
save(gs_gencode_v25_hg38, file = 'gs_gencode_v25_hg38.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
