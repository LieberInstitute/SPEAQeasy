FROM libddocker/bcftools-1.3.1:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
RUN apt-get install -y libncurses5-dev

## Define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

### vcftools
WORKDIR $SRC
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz && \
    tar -xjvf vcftools-0.1.15.tar.gz && \
    cd vcftools-0.1.15 && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    cp -r * $BIN