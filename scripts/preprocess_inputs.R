library('BiocParallel')
library('devtools')
library('jaffelab')

run_command = function(command) {
    print(paste0("Running command: '", command, "'..."))
    system(command)
    print("Done.")
}

manifest <- read.table("samples.manifest", header = FALSE, stringsAsFactors = FALSE)

## Is the data paired end?
paired <- ncol(manifest) > 3

############################################################
#  Verify the files in the manifest have known extensions
############################################################

print("Verifying file extensions from the manifest are valid...")

#  Get the FASTQ filenames, as declared in samples.manifest
filenames = manifest[,1]
if (paired) filenames = c(filenames, manifest[,3])

is_zipped = grepl(".gz", filenames, fixed=TRUE)

valid_exts = c('fastq.gz', 'fq.gz', 'fastq', 'fq')
actual_exts = sapply(1:length(filenames), function(i) {
    temp = strsplit(sub(".gz", "", filenames[i]), ".", fixed=TRUE)[[1]]
    if (is_zipped[i]) {
        return(paste0(temp[length(temp)], '.gz'))
    } else {
        return(temp[length(temp)])
    }
})

if (!all(actual_exts %in% valid_exts)) {
    stop("Unrecognized fastq filename extension. Should be fastq.gz, fq.gz, fastq or fq")
}

if (paired && any(actual_exts[1:nrow(manifest)] != actual_exts[(nrow(manifest)+1):length(actual_exts)])) {
    stop("A given pair of reads must have the same file extensions.")
}

####################################################################
#  Perform merging and renaming of files for use in the
#  pipeline. Then write a new manifest to reflect these actions
####################################################################

#  This forms a list, where each element is a vector containing row numbers
#  of the manifest to combine (and each element contains a unique set of rows)
indicesToCombine = list()
for (i in 1:nrow(manifest)) {
    idMatchRows = which(manifest[, ncol(manifest)] == manifest[i, ncol(manifest)])
    if (all(idMatchRows >= i) && length(idMatchRows) > 1) {
        indicesToCombine[[i]] = idMatchRows
    }
}
indicesToCombine = indicesToCombine[!sapply(indicesToCombine, is.null)]

if (paired) {
    #  Merge files that require merging
    print("Merging any files that need to be merged...")
    for (indices in indicesToCombine) {
        #  Determine file extension for the merged file to have
        first_ext = actual_exts[indices[1]]
        
        #  Do the file merging
        files_to_combine = do.call(paste, as.list(manifest[indices, 1]))
        new_file = paste0(manifest[indices[1], 5], '_1.', first_ext)
        command = paste('cat', files_to_combine, '>', new_file)
        run_command(command)
        
        files_to_combine = do.call(paste, as.list(manifest[indices, 3]))
        new_file = paste0(manifest[indices[1], 5], '_2.', first_ext)
        command = paste('cat', files_to_combine, '>', new_file)
        run_command(command)
    }
    
    #  Symbolically link any remaining files: this renames the files
    #  by their associated sampleID, and uses the paired suffices _1 and _2.
    print("Renaming and symbolically linking files for handling in the pipeline...")
    if (length(unlist(indicesToCombine)) > 0) {
        remaining_rows = (1:nrow(manifest))[-unlist(indicesToCombine)]
    } else {
        remaining_rows = 1:nrow(manifest)
    }
    for (index in remaining_rows) {
        first_ext = actual_exts[index]
        
        new_file = paste0(manifest[index, 5], '_1.', first_ext)
        command = paste('ln -s', manifest[index, 1], new_file)
        run_command(command)
        
        new_file = paste0(manifest[index, 5], '_2.', first_ext)
        command = paste('ln -s', manifest[index, 3], new_file)
        run_command(command)
    }
    
    #  Rewrite a manifest to reflect the file name changes and any merging
    print("Constructing a new manifest to reflect these changes...")
    first_reads = basename(system('ls *_1.f*q*', intern=TRUE))
    ids = ss(first_reads, '_1.', fixed=TRUE)
    new_man = paste(first_reads,
                    0,
                    basename(system('ls *_2.f*q*', intern=TRUE)),
                    0,
                    ids)
    writeLines(new_man, con="samples_processed.manifest")

} else {
    print("Merging any files that need to be merged...")
    for (indices in indicesToCombine) {
        #  Do the file merging
        files_to_combine = do.call(paste, as.list(manifest[indices, 1]))
        new_file = paste0(manifest[indices[1], 3], '.', actual_exts[indices[1]])
        command = paste('cat', files_to_combine, '>', new_file)
        run_command(command)
    }
    
    #  Symbolically link any remaining files; also renames the files
    #  by their associated sampleID
    print("Renaming and symbolically linking files for handling in the pipeline...")
    if (length(unlist(indicesToCombine)) > 0) {
        remaining_rows = (1:nrow(manifest))[-unlist(indicesToCombine)]
    } else {
        remaining_rows = 1:nrow(manifest)
    }
    for (index in remaining_rows) {
        new_file = paste0(manifest[index, 3], '.', actual_exts[index])
        command = paste('ln -s', manifest[index, 1], new_file)
        run_command(command)
    }
    
    #  Rewrite a manifest to reflect the file name changes and any merging
    print("Constructing a new manifest to reflect these changes...")
    reads = basename(system('ls *.f*q*', intern=TRUE))
    ids = ss(reads, '.', fixed=TRUE)
    new_man = paste(reads, 0, ids)
    writeLines(new_man, con="samples_processed.manifest")
}

print("Done all tasks.")
