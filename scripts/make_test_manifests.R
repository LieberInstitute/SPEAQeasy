library('jaffelab')

for (species in c("human", "mouse", "rat")) {
    for (pairing in c("single", "paired")) {
        for (strand in c("stranded", "unstranded")) {
            man_dir = file.path(getwd(), "test", species, pairing, strand)
            
            if (pairing == "single") {
                fq_files = file.path(man_dir, list.files(man_dir, ".*\\.fastq($|\\.gz$)"))
                sample_names = basename(ss(gsub("_file[12]", "", fq_files), "\\.fastq"))
                
                man_lines = paste(fq_files, 0, sample_names, sep="\t")
            } else {
                fq_files1 = file.path(man_dir, list.files(man_dir, ".*_1\\.fastq($|\\.gz$)"))
                fq_files2 = file.path(man_dir, list.files(man_dir, ".*_2\\.fastq($|\\.gz$)"))
                sample_names = basename(ss(gsub("(_file[12])*_1", "", fq_files1), "\\.fastq"))
                
                man_lines = paste(fq_files1, 0, fq_files2, 0, sample_names, sep="\t")
            }
            
            writeLines(man_lines, paste0(man_dir, "/samples.manifest"))
        }
    }
}
