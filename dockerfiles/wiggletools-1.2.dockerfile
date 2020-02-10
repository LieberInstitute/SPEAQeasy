FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
RUN apt-get install bc

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

### LibBigWig
WORKDIR $SRC
RUN git clone https://github.com/dpryan79/libBigWig.git && \
    cd libBigWig && \
    make install && \
    cp -r * $BIN

### WiggleTools
WORKDIR $SRC
RUN git clone https://github.com/Ensembl/WiggleTools.git && \
    cd WiggleTools && \
    make && \
    cd bin && \
    cp -r * $BIN

### wigToBigWig install 
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    chmod 755 wigToBigWig && \
    cp wigToBigWig $BIN

WORKDIR $BIN
ENV LD_LIBRARY_PATH /usr/local/lib
