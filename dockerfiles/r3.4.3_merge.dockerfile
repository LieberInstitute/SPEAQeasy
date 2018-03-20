FROM libd/r3.4.3_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update Core ########
RUN apt-get update

######## Install R packages:versions ########
## getopt:1.20.1
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/getopt_1.20.1.tar.gz',type='source',dependencies=TRUE)"

## devtools:1.13.4
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/devtools_1.13.4.tar.gz',type='source',dependencies=TRUE)"
 
## BiocParallel:1.12.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/BiocParallel_1.12.0.tar.gz',type='source',dependencies=TRUE)"

