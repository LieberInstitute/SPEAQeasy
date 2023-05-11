#  Given two 'count_objects' directories from two SPEAQeasy runs on identical
#  samples, verify that all objects are identical between the runs (excluding
#  variables that refer to full paths)

#  We assume there are 6 outputs directories in total (two sets of runs for
#  each reference genome, where the second set typically occurs after some
#  significant update to SPEAQeasy's code), structured and named like this:
#  [base_dir]/out_["old" or "new"]_[reference_genome]

library("SummarizedExperiment")
library("getopt")

spec <- matrix(
    c(
        "base_dir", "d", 1, "character", "Directory containing outputs for each run",
        "stop_on_err", "e", 1, "logical", "Halt upon any difference in outputs?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

for (ref in c("hg38", "mm10", "rn6")) {
    out_dir1 <- file.path(opt$base_dir, paste0("out_old_", ref), "count_objects")
    out_dir2 <- file.path(opt$base_dir, paste0("out_new_", ref), "count_objects")

    rdas1 <- list.files(out_dir1, pattern = ".*\\.(rda|Rdata)$", full.names = TRUE)
    rdas2 <- list.files(out_dir2, pattern = ".*\\.(rda|Rdata)$", full.names = TRUE)

    print(paste0("Comparing Rdata/rda files between sets for ", ref, "..."))
    for (i in 1:length(rdas1)) {
        current_objects <- ls()
        load(rdas1[i])
        diff_objects <- ls()[!(ls() %in% c(current_objects, "current_objects"))]

        obj_list1 <- list()
        for (diff_obj in diff_objects) {
            obj_list1[[diff_obj]] <- eval(parse(text = diff_obj))
        }

        load(rdas2[i])

        obj_list2 <- list()
        for (diff_obj in diff_objects) {
            obj_list2[[diff_obj]] <- eval(parse(text = diff_obj))
            rm(list = diff_obj)
        }

        for (j in 1:length(obj_list1)) {
            if (!identical(obj_list1[[j]], obj_list2[[j]])) {
                #  If a data frame, allow different odering of columns
                if (class(obj_list1[[j]]) == "data.frame") {
                    stopifnot(all(colnames(obj_list1[[j]]) %in% colnames(obj_list2[[j]])))
                    stopifnot(all(colnames(obj_list2[[j]]) %in% colnames(obj_list1[[j]])))

                    obj_list1[[j]] <- obj_list1[[j]][, colnames(obj_list2[[j]])]

                    if (identical(obj_list1[[j]], obj_list1[[j]])) {
                        print("Found a pair of data frames that were identical only after re-ordering columns.")
                    } else {
                        out_str <- paste0("Not all objects are identical between files named '", basename(rdas1[i]), "'.")
                        if (opt$stop_on_err) {
                            stop(out_str)
                        } else {
                            warning(out_str)
                        }
                    }
                } else if (class(obj_list1[[j]]) == "RangedSummarizedExperiment") {
                    if (names(obj_list1[j]) %in% c("rse_exon", "rse_gene", "rse_tx", "rse_jx")) {
                        #  Remove bamFile column in colData
                        colData(obj_list1[[j]])$bamFile <- NULL
                        colData(obj_list2[[j]])$bamFile <- NULL

                        #  Allow re-ordering of columns in rowData
                        stopifnot(all(colnames(rowData(obj_list1[[j]])) %in% colnames(rowData(obj_list2[[j]]))))
                        stopifnot(all(colnames(rowData(obj_list2[[j]])) %in% colnames(rowData(obj_list1[[j]]))))

                        rowData(obj_list1[[j]]) <- rowData(obj_list1[[j]])[, colnames(rowData(obj_list2[[j]]))]
                        
                        #  Assume that we're comparing outputs for the sake of
                        #  testing, and that during testing SPEAQeasy might be
                        #  run with slightly different settings (which we don't
                        #  care about)
                        metadata(obj_list1[[j]])$SPEAQeasy_settings = NULL
                        metadata(obj_list2[[j]])$SPEAQeasy_settings = NULL
                    }
                    if (identical(obj_list1[[j]], obj_list2[[j]])) {
                        print("Removed the 'bamFile' columns and SPEAQeasy settings from two RSEs, after which point the objects were equivalent.")
                    } else {
                        out_str <- paste0("Not all objects are identical between files named '", basename(rdas1[i]), "'.")
                        if (opt$stop_on_err) {
                            stop(out_str)
                        } else {
                            warning(out_str)
                        }
                    }
                }
            }
        }
    }

    #  Compare CSV outputs, excluding the 'bamFile' column which should have
    #  different values
    csv1 <- read.csv(list.files(out_dir1, pattern = ".*\\.csv$", full.names = TRUE))
    csv1$bamFile <- NULL

    csv2 <- read.csv(list.files(out_dir2, pattern = ".*\\.csv$", full.names = TRUE))
    csv2$bamFile <- NULL

    if (!identical(csv1, csv2)) {
        stop("Metrics CSVs were not identical.")
    }

    if (opt$stop_on_err) {
        print(paste0("All objects were identical for ", ref, "."))
    }
}

print("All objects across all references were identical!")
