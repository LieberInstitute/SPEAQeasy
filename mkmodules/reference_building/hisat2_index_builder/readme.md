# mkmodule: hisat2_index_builder

## Description

This mkmodule downloads a .fa.gz file, uncompresses it, and uses it to build a hisat2 index.

## Expected outputs

1. Many .ht2 hisat2 index files for the given reference
2. .fa.gz file for the given reference genome
3. .fa file for the given reference genome

## Software Dependencies

For this module to work, make sure your $PATH includes access to the following executables:

- mk (from https://9fans.github.io/plan9port/man/man1/install.html) see man here: https://9fans.github.io/plan9port/man/man1/mk.html
- hisat2-build (from https://ccb.jhu.edu/software/hisat2/manual.shtml#building-from-source)
- gunzip
- wget

## Variables Required

The following variables should be specified when invoking this module (via bash runmk.sh VARIABLE1="MY_VALUE", VARIABLE2="MY_OTHER_VALUE, etc.")

* REFERENCE_FASTA_URL <- full URL for compressed fasta download

* NUMBER_OF_THREADS <- number of threads that will be used during hisat2 index building

**NOTE**: if this values are not provided as indicated in the execution Example, default values will be used taken from the config.mk file in this module. config.mk file can be used for quick testing and developing.

## Execution Example:
` bash runmk.sh REFERENCE_FASTA_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/CHR_22/hs_ref_GRCh38.p12_chr22.fa.gz" NUMBER_OF_THREADS="4" `

## Author information
Developed by Israel Aguilar (iaguilaror@gmail.com) on NOV 2018
