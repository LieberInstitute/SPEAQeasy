library('jaffelab')

man = read.table('samples_processed.manifest', header=FALSE, stringsAsFactors=FALSE)

strand_files = list.files(pattern='.*_strandness_pattern\\.txt')
strands = sapply(strand_files, readLines)
ids = ss(strand_files, '_strandness', fixed=TRUE)

strand_col = strands[match(ids, man[,ncol(man)])]
man = cbind(man, strand_col)

write.table(man, 'samples_complete.manifest', row.names=FALSE, col.names=FALSE,
            quote=FALSE)
