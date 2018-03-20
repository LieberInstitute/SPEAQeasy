FROM libddocker/ubuntu16.04_base:latest
 
MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"
 
########### Update Core ###########
RUN apt-get update
 
## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install cmake
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:george-edison55/cmake-3.x
RUN apt-get update
RUN apt-get install -y cmake
RUN apt-get upgrade -y

## Install zlib, HDF5 libraries
RUN apt-get install -y \
    zlib1g-dev \
    libhdf5-dev \
    autoconf \
    autotools-dev

## Install Kallisto
WORKDIR $SRC
RUN git clone https://github.com/pachterlab/kallisto.git && \
    cd kallisto && \
    cd ext/htslib && \
    autoheader && \
    autoconf && \
    cd ../.. && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    cp -r * $BIN
