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
library("getopt")

spec <- matrix(c("repo_dir", "d", 1, "character", "path to SPEAQeasy repo"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

for (species in c("human", "mouse", "rat")) {
    for (pairing in c("single", "paired")) {
        for (strand in c("unstranded", "forward", "reverse")) {
            man_dir <- file.path("test", species, pairing, strand)

            if (pairing == "single") {
                fq_files <- file.path(opt$repo_dir, man_dir, list.files(man_dir, ".*\\.fastq($|\\.gz$)"))
                sample_names <- basename(ss(gsub("_file[12]", "", fq_files), "\\.fastq"))

                man_lines <- paste(fq_files, 0, sample_names, sep = "\t")
            } else {
                fq_files1 <- file.path(opt$repo_dir, man_dir, list.files(man_dir, ".*_1\\.fastq($|\\.gz$)"))
                fq_files2 <- file.path(opt$repo_dir, man_dir, list.files(man_dir, ".*_2\\.fastq($|\\.gz$)"))
                sample_names <- basename(ss(gsub("(_file[12])*_1", "", fq_files1), "\\.fastq"))

                man_lines <- paste(fq_files1, 0, fq_files2, 0, sample_names, sep = "\t")
            }

            writeLines(man_lines, paste0(man_dir, "/samples.manifest"))
        }
    }
}
