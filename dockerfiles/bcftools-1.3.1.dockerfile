FROM libddocker/samtools-1.3.1:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install bcftools
WORKDIR $SRC
RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.3.1 && \
    make && \
    make prefix=/usr/local/bin install && \
    cp -r * $BIN