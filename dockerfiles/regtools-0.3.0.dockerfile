FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update 

## Install cmake
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:george-edison55/cmake-3.x
RUN apt-get update
RUN apt-get install -y cmake
RUN apt-get upgrade -y

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
