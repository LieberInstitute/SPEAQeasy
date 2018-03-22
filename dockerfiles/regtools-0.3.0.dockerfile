FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update 

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install Regtools
WORKDIR $SRC
RUN git clone https://github.com/griffithlab/regtools && \
    cd regtools && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cp -r * $BIN
