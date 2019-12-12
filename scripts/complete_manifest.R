library('jaffelab')

man = read.table('samples_processed.manifest', header=FALSE, stringsAsFactors=FALSE)

strand_files = list.files('.*\.strandness_pattern\.txt')
strands = sapply(strand_files, readLines)
ids = ss(strand_files, '.strandness', fixed=TRUE)

strand_col = strands[match(ids, man[,ncol(man)])]
man = cbind(man, strand_col)

write.table(man, 'samples_complete.manifest', sep=' ', header=FALSE, stringsAsFactors=FALSE)