FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

############ Update Core ############
RUN apt-get update

## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install FastQC
WORKDIR $SRC
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    chmod -R 755 /usr/local/src/FastQC && \
    ln -s /usr/local/src/FastQC/fastqc $BIN
