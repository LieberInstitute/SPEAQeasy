packages <- c('Biostrings', 'GenomicRanges', 'GenomicFeatures', 'org.Hs.eg.db',
    'biomaRt', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
    'org.Mm.eg.db', 'BSgenome.Mmusculus.UCSC.mm10', 'org.Rn.eg.db',
    'BSgenome.Rnorvegicus.UCSC.rn6', 'derfinder', 'bumphunter', 'jaffelab',
    'devtools', 'getopt', 'BiocParallel', 'rafalib', 'SummarizedExperiment',
    'plyr', 'rtracklayer', 'RColorBrewer')

## Try to load them
load_res <- sapply(packages, requireNamespace, quietly = TRUE)

if(any(!load_res)) {
    system('touch .missing_R_packages')
}
