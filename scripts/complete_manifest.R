library("here")

#  If R was installed locally, use packages associated with that installation
if (.libPaths()[1] == here("Software", "R-4.1.0", "library")) {
    library("checkpoint")
    checkpoint("2021-09-01",
        project_dir = here("scripts", "r_packages"),
        checkpoint_location = here("Software", "R-4.1.0")
    )
}

library("jaffelab")

man <- read.table("samples_processed.manifest", header = FALSE, stringsAsFactors = FALSE)

strand_files <- list.files(pattern = ".*_strandness_pattern\\.txt")
strands <- sapply(strand_files, readLines)
ids <- ss(strand_files, "_strandness", fixed = TRUE)

strand_col <- strands[match(ids, man[, ncol(man)])]
man <- cbind(man, strand_col)

write.table(man, "samples_complete.manifest",
    row.names = FALSE, col.names = FALSE,
    quote = FALSE
)
