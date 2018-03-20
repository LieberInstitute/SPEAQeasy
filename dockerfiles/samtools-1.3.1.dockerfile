FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
RUN apt-get install -y libncurses5-dev

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install Samtools
WORKDIR $SRC
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-1.3.1 && \
    make && \
    make prefix=/usr/local/bin install && \
    cp -r * $BIN
