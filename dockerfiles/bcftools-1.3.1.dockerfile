FROM libddocker/samtools-1.3.1:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
RUN apt-get install -y libncurses5-dev

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

### htslib
WORKDIR $SRC
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && \
    make install && \
    cp -r * $BIN

## Install bcftools
WORKDIR $SRC
RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.3.1 && \
    make && \
    make prefix=/usr/local/bin install && \
    cp -r * $BIN