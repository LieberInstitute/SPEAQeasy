library("Biostrings")
library("jaffelab")

chrom_names <- paste0(c(1:20, "X", "Y", "MT"))
seq_width_chars <- 80

fasta <- readDNAStringSet(list.files(pattern = ".*\\.fa$"))

#   Subset to canonical reference chromosomes
chr <- ss(ss(names(fasta), " ", 3), ":", 3)
fasta <- fasta[chr %in% chrom_names]

#   Reconstruct the lines of the output FASTA from the sequence headers and
#   actual sequence
out_lines <- c()
# for (i in 1:length(fasta)) {
#     #   Add sequence header
#     out_lines = c(out_lines, paste0('>', names(fasta)[i]))
#
#     #   Split the sequence into chunks of [seq_width_chars] per line
#     this_line = strsplit(
#         gsub(
#             paste0("([ACTGN]{", seq_width_chars, "})"),
#             "\\1 ",
#             as.character(fasta[[i]])
#         ),
#         " "
#     )[[1]]
#
#     #   Add this sequence to the FASTA lines
#     out_lines = c(out_lines, this_line)
# }

#   For computational speed, print the entire sequence on one line (breaking
#   the convention to use 80 bps per line)
for (i in 1:length(fasta)) {
    out_lines <- c(
        out_lines,
        paste0(">", names(fasta)[i]),
        as.character(fasta[[i]])
    )
}

#   Overwrite the full FASTA
writeLines(out_lines, con = list.files(pattern = ".*\\.fa$"))
