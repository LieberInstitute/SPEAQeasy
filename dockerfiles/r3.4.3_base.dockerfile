FROM r-base:3.4.3

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ########
RUN apt-get update --fix-missing
RUN apt-get install -y -f  apt-utils
RUN apt-get install -y -f libcurl4-gnutls-dev
RUN apt-get install -y -f libcurl4-openssl-dev
RUN apt-get install -y -f libexpat1-dev
RUN apt-get install -y -f libpython-dev
RUN apt-get update --fix-missing
RUN apt-get install -y -f libglib2.0-bin
RUN apt-get install -y -f libglib2.0-dev
RUN apt-get install -y -f libcairo2=1.15.10-1 --allow-downgrades
RUN apt-get install -y \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libxml2-dev \
    default-jdk \
    mesa-common-dev \
    libglu1-mesa \
    libxi-dev \
    libxmu-dev \
    vim \
    wget \
    curl \
    git \
    zip \
    sudo
RUN apt-get clean
RUN apt-get update

######## Install R Base Packages ########
RUN apt-get update
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/BiocInstaller_1.28.0.tar.gz',type='source',dependencies=TRUE)"

RUN Rscript -e "install.packages(c('remotes','devtools','RColorBrewer','plyr'),repo='http://cran.rstudio.com/',dependencies=TRUE)"

## getopt:1.20.2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/getopt_1.20.2.tar.gz',type='source',dependencies=TRUE)"

## RCurl: 1.95-4.10
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/RCurl_1.95-4.10.tar.gz',type='source',dependencies=TRUE)"

## XML: 3.98-1.10
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/XML_3.98-1.10.tar.gz',type='source',dependencies=TRUE)"

## zlibbioc: 1.24.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/zlibbioc_1.24.0.tar.gz',type='source',dependencies=TRUE)"

## BiocGenerics: 0.24.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/BiocGenerics_0.24.0.tar.gz',type='source',dependencies=TRUE)"

## S4Vectors: 0.16.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/S4Vectors_0.16.0.tar.gz',type='source',dependencies=TRUE)"

## IRanges: 2.12.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/IRanges_2.12.0.tar.gz',type='source',dependencies=TRUE)"

## XVector: 0.18.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/XVector_0.18.0.tar.gz',type='source',dependencies=TRUE)"

## Biostrings: 2.46.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/Biostrings_2.46.0.tar.gz',type='source',dependencies=TRUE)"

## GenomeInfoDbData: 1.0.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/data/annotation/src/contrib/GenomeInfoDbData_1.0.0.tar.gz',type='source',dependencies=TRUE)"

## GenomeInfoDb: 1.14.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/GenomeInfoDb_1.14.0.tar.gz',type='source',dependencies=TRUE)"

## Genomic Ranges: 1.30.3
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/GenomicRanges_1.30.3.tar.gz',type='source',dependencies=TRUE)"

## lambda.r: 1.2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/lambda.r_1.2.tar.gz',type='source',dependencies=TRUE)"

## futile.options: 1.0.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/futile.options_1.0.0.tar.gz',type='source',dependencies=TRUE)"

## futile.logger: 1.4.3
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/futile.logger_1.4.3.tar.gz',type='source',dependencies=TRUE)"

## snow: 0.4-2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/snow_0.4-2.tar.gz',type='source',dependencies=TRUE)"

## BH: 1.66.0-1
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/BH_1.66.0-1.tar.gz',type='source',dependencies=TRUE)"

## BiocParallel: 1.12.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/BiocParallel_1.12.0.tar.gz',type='source',dependencies=TRUE)"

## Rsamtools: 1.30.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/Rsamtools_1.30.0.tar.gz',type='source',dependencies=TRUE)"

## BioBase: 2.38.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/Biobase_2.38.0.tar.gz',type='source',dependencies=TRUE)"

## maxtrixStats: 0.53.1
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/matrixStats_0.53.1.tar.gz',type='source',dependencies=TRUE)"

## DelayedArray: 0.4.1
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/DelayedArray_0.4.1.tar.gz',type='source',dependencies=TRUE)"

## Summarize Experiment: 1.8.1 
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/SummarizedExperiment_1.8.1.tar.gz',type='source',dependencies=TRUE)"

## GenomicAlignments: 1.14.1
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/GenomicAlignments_1.14.1.tar.gz',type='source',dependencies=TRUE)"

## rtracklayer: 1.38.3
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/rtracklayer_1.38.3.tar.gz',type='source',dependencies=TRUE)"

## Biobase: 2.38.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/Biobase_2.38.0.tar.gz',type='source',dependencies=TRUE)"

## DBI: 0.8
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/DBI_0.8.tar.gz',type='source',dependencies=TRUE)"

## bit: 1.1-12
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/bit_1.1-12.tar.gz',type='source',dependencies=TRUE)"

## bit64: 0.9.7
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/bit64_0.9-7.tar.gz',type='source',dependencies=TRUE)"

## crayon: 1.3.4
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/crayon_1.3.4.tar.gz',type='source',dependencies=TRUE)"

## assertthat: 0.2.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/assertthat_0.2.0.tar.gz',type='source',dependencies=TRUE)"

## utf8: 1.1.3
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/utf8_1.1.3.tar.gz',type='source',dependencies=TRUE)"

## rlang: 0.2.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/rlang_0.2.0.tar.gz',type='source',dependencies=TRUE)"

## cli: 1.0.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/cli_1.0.0.tar.gz',type='source',dependencies=TRUE)"

## pillar: 1.2.1
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/pillar_1.2.1.tar.gz',type='source',dependencies=TRUE)"

## tibble: 1.4.2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/tibble_1.4.2.tar.gz',type='source',dependencies=TRUE)"

## blob: 1.1.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/blob_1.1.0.tar.gz',type='source',dependencies=TRUE)"

## plogr: 0.1-1
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/plogr_0.1-1.tar.gz',type='source',dependencies=TRUE)"

## RSQLite: 2.0
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/RSQLite_2.0.tar.gz',type='source',dependencies=TRUE)"

## AnnotationDbi: 1.40.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/AnnotationDbi_1.40.0.tar.gz',type='source',dependencies=TRUE)"

## foreach: 1.4.4
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/foreach_1.4.4.tar.gz',type='source',dependencies=TRUE)"

## iterators: 1.0.9
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/iterators_1.0.9.tar.gz',type='source',dependencies=TRUE)"

## locfit: 1.5-9.1
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/locfit_1.5-9.1.tar.gz',type='source',dependencies=TRUE)"

## limma: 3.34.9
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/limma_3.34.9.tar.gz',type='source',dependencies=TRUE)"

## registry: 0.5
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/registry_0.5.tar.gz',type='source',dependencies=TRUE)"

## xtable: 1.8-2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/xtable_1.8-2.tar.gz',type='source',dependencies=TRUE)"

## pkgmaker: 0.22
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/pkgmaker_0.22.tar.gz',type='source',dependencies=TRUE)"

## rngtools: 1.2.4
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/rngtools_1.2.4.tar.gz',type='source',dependencies=TRUE)"

## doRNG: 1.6.6
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/doRNG_1.6.6.tar.gz',type='source',dependencies=TRUE)"

## knitr: 1.20
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/knitr_1.20.tar.gz',type='source',dependencies=TRUE)"

## testthat: 2.0.0
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/testthat_2.0.0.tar.gz',type='source',dependencies=TRUE)"

## magrittr: 1.5
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/magrittr_1.5.tar.gz',type='source',dependencies=TRUE)"

## glue: 1.2.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/glue_1.2.0.tar.gz',type='source',dependencies=TRUE)"

## stringr: 1.3.0
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/stringr_1.3.0.tar.gz',type='source',dependencies=TRUE)"

## prettyunits: 1.0.2
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/prettyunits_1.0.2.tar.gz',type='source',dependencies=TRUE)"

## R6
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/R6_2.2.2.tar.gz',type='source',dependencies=TRUE)"

## progress
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/progress_1.1.2.tar.gz',type='source',dependencies=TRUE)"

## bitops: 1.0-6
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/bitops_1.0-6.tar.gz',type='source',dependencies=TRUE)"

## XML: 3.98-1.10
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/XML_3.98-1.10.tar.gz',type='source',dependencies=TRUE)"

## Rcurl: 1.95-4.10
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/RCurl_1.95-4.10.tar.gz',type='source',dependencies=TRUE)"

## derfinderHelper: 1.12.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/derfinderHelper_1.12.0.tar.gz',type='source',dependencies=TRUE)"

## install mysql client
RUN apt-get install -y libmariadb-client-lgpl-dev
RUN apt-get update

## RMySQL: 0.10.14
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/RMySQL_0.10.14.tar.gz',type='source',dependencies=TRUE)"

## biomaRt: 2.34.2
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/biomaRt_2.34.2.tar.gz',type='source',dependencies=TRUE)"

## GenomicFeatures: 1.30.3
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/GenomicFeatures_1.30.3.tar.gz',type='source',dependencies=TRUE)"

## bumphunter: 1.20.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/bumphunter_1.20.0.tar.gz',type='source',dependencies=TRUE)"

## zlibbioc: 1.24.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/zlibbioc_1.24.0.tar.gz',type='source',dependencies=TRUE)"

## Rsamtools: 1.30.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/Rsamtools_1.30.0.tar.gz',type='source',dependencies=TRUE)"

## BSgenome: 1.46.0
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/BSgenome_1.46.0.tar.gz',type='source',dependencies=TRUE)"

## VariantAnnotation: 1.24.5
RUN Rscript -e "install.packages('https://bioconductor.org/packages/release/bioc/src/contrib/VariantAnnotation_1.24.5.tar.gz',type='source',dependencies=TRUE)"

## GenomicFiles: 1.14.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/GenomicFiles_1.14.0.tar.gz',type='source',dependencies=TRUE)"

## Formula: 1.2-2
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/Formula_1.2-2.tar.gz',type='source',dependencies=TRUE)"

## survival: 2.41-3
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/survival_2.41-3.tar.gz',type='source',dependencies=TRUE)"

## gtable: 0.2.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/gtable_0.2.0.tar.gz',type='source',dependencies=TRUE)"

## dichromatL 2.0-0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/dichromat_2.0-0.tar.gz',type='source',dependencies=TRUE)"

## labeling: 0.3.
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/labeling_0.3.tar.gz',type='source',dependencies=TRUE)"

## colorspace
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/colorspace_1.3-2.tar.gz',type='source',dependencies=TRUE)"

## munsell: 0.4.3
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/munsell_0.4.3.tar.gz',type='source',dependencies=TRUE)"

## viridisLite: 0.3.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/viridisLite_0.3.0.tar.gz',type='source',dependencies=TRUE)"

## scales: 0.5.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/scales_0.5.0.tar.gz',type='source',dependencies=TRUE)"

## reshape2: 1.4.3
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/reshape2_1.4.3.tar.gz',type='source',dependencies=TRUE)"

## ggplot2: 2.2.1
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/ggplot2_2.2.1.tar.gz',type='source',dependencies=TRUE)"

## lattice: 0.20-35
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/lattice_0.20-35.tar.gz',type='source',dependencies=TRUE)"

## RColorBrewer: 1.1-2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/RColorBrewer_1.1-2.tar.gz',type='source',dependencies=TRUE)"

## latticeExtra: 0.6-28
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/latticeExtra_0.6-28.tar.gz',type='source',dependencies=TRUE)"

## acepack: 1.4.1
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/acepack_1.4.1.tar.gz',type='source',dependencies=TRUE)"

## gridExtra: 2.3
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/gridExtra_2.3.tar.gz',type='source',dependencies=TRUE)"

## viridis: 0.5.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/viridis_0.5.0.tar.gz',type='source',dependencies=TRUE)"

## checkmate: 1.8.5
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/checkmate_1.8.5.tar.gz',type='source',dependencies=TRUE)"

## htmlwidgets: 1.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/htmlwidgets_1.0.tar.gz',type='source',dependencies=TRUE)"

## htmlTable: 1.11.2
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/htmlTable_1.11.2.tar.gz',type='source',dependencies=TRUE)"

## data.table: 1.10.4-3
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/data.table_1.10.4-3.tar.gz',type='source',dependencies=TRUE)"

## Hmisc: 4.1-1
RUN Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/Hmisc_4.1-1.tar.gz',type='source',dependencies=TRUE)"

## qvalue: 2.10.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/qvalue_2.10.0.tar.gz',type='source',dependencies=TRUE)"

## derfinder: 1.12.6
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/bioc/src/contrib/derfinder_1.12.6.tar.gz',type='source',dependencies=TRUE)"

## jaffelab package
RUN Rscript -e "devtools::install_github('LieberInstitute/jaffelab')" 

## org.Hs.eg.db: 3.5.0
RUN Rscript -e "install.packages('http://bioconductor.org/packages/release/data/annotation/src/contrib/org.Hs.eg.db_3.5.0.tar.gz',type='source',dependencies=TRUE)"

## rafalib: 1.0.0
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/rafalib_1.0.0.tar.gz',type='source',dependencies=TRUE)"

## plyr: 1.8.4
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/plyr_1.8.4.tar.gz',type='source',dependencies=TRUE)"