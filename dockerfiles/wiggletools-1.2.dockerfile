FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
RUN apt-get install -y \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libboost-all-dev \
    libncurses5-dev \
    python-dev \
    libbz2-dev

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

### htslib
WORKDIR $SRC
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-{version} && \
    make && \
    make install && \
    cp -r * $BIN

### LibBigWig
WORKDIR $SRC
RUN git clone https://github.com/dpryan79/libBigWig.git && \
    cd libBigWig && \
    make install && \
    cp -r * $BIN

### GSL Lib
WORKDIR $SRC
RUN wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && \
    tar -xvzpf gsl-latest.tar.gz && \
    cd gsl* && \
    ./configure && \
    make && \
    make install && \
    cp -r * $BIN

### WiggleTools
WORKDIR $SRC
RUN git clone https://github.com/Ensembl/WiggleTools.git && \
    cd WiggleTools && \
    make && \
    cd bin && \
    cp -r * $BIN

### wigToBigWig install 
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    chmod 755 wigToBigWig && \
    cp wigToBigWig $BIN

WORKDIR $BIN
ENV LD_LIBRARY_PATH /usr/local/lib
