FROM libddocker/ubuntu16.04_base:latest
 
MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"
 
########### Update Core ###########
RUN apt-get update
 
## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

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
