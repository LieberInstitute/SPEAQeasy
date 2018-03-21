FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install SubRead
WORKDIR $SRC
RUN wget https://sourceforge.net/projects/subread/files/subread-1.6.0/subread-1.6.0-source.tar.gz && \
    tar -zxvf subread-1.6.0-source.tar.gz && \
    cd subread-1.6.0-source/src && \
    make -f Makefile.Linux && \
    cd ../bin && \
    cp -r * $BIN
WORKDIR $BIN

