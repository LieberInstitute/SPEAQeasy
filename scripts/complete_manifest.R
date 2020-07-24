#  If R was installed locally, use packages associated with that installation
if (grepl("SPEAQeasy/Software/R-3.6.1/library", .libPaths()[1])) {
    repo_dir = dirname(dirname(dirname(.libPaths()[1])))
    library('checkpoint')
    checkpoint::checkpoint("2019-08-05",
                           project=paste0(repo_dir, '/scripts/'),
                           checkpointLocation = paste0(repo_dir, '/Software/R-3.6.1/'))
}

library('jaffelab')

man = read.table('samples_processed.manifest', header=FALSE, stringsAsFactors=FALSE)

strand_files = list.files(pattern='.*_strandness_pattern\\.txt')
strands = sapply(strand_files, readLines)
ids = ss(strand_files, '_strandness', fixed=TRUE)

strand_col = strands[match(ids, man[,ncol(man)])]
man = cbind(man, strand_col)

write.table(man, 'samples_complete.manifest', row.names=FALSE, col.names=FALSE,
            quote=FALSE)
