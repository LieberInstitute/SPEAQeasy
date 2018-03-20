FROM libddocker/ubuntu16.04_base:latest
 
MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"
 
########### Update Core ###########
RUN apt-get update
 
## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install HISAT2
WORKDIR $SRC
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip && \ 
    unzip hisat2-2.0.4-Linux_x86_64.zip && \
    chmod -R 755 hisat2-2.0.4 && \
    cp -r hisat2-2.0.4/* $BIN

WORKDIR $BIN
