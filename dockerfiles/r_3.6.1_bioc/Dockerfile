FROM bioconductor/bioconductor_full:RELEASE_3_10

COPY check_R_packages.R ./check_R_packages.R
RUN Rscript check_R_packages.R
