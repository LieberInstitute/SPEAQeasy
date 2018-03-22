FROM ubuntu:16.04

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

########### Update Core ###########
RUN apt-get update
 
## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

######## Install Core ######## 
RUN apt-get update
RUN apt-get install -y software-properties-common \
    python-software-properties \
    libreadline6 \
    libreadline6-dev
RUN add-apt-repository ppa:george-edison55/cmake-3.x
RUN apt-get update --fix-missing
RUN apt-get install -y -f apt-utils
RUN apt-get install -y -f libcurl4-gnutls-dev
RUN apt-get install -y -f libcurl4-openssl-dev
RUN apt-get install -y -f libexpat1-dev
RUN apt-get install -y -f libpython-dev
RUN apt-get update --fix-missing
RUN apt-get install -y \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libxml2-dev \
    default-jdk \
    mesa-common-dev \
    libglu1-mesa \
    libxi-dev \
    libxmu-dev \
    vim \
    wget \
    curl \
    git \
    zip \
    sudo \
    python-pip \
    python-setuptools \
    python-dev \
    bzip2 \
    zlib1g-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libboost-all-dev \
    libncurses5-dev \
    python-dev \
    libbz2-dev \
    libncurses5-dev \
    autoconf \
    autogen \
    build-essential \
    cmake-curses-gui \
    liblzo2-dev \
    zlib1g-dev \
    libhdf5-dev \
    autoconf \
    autotools-dev
RUN apt-get clean


### htslib
WORKDIR $SRC
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.3.2 && \
    make && \
    make install && \
    cp -r * $BIN

### GSL Lib
WORKDIR $SRC
RUN wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && \
    tar -zxvf gsl-latest.tar.gz && \
    cd gsl* && \
    ./configure && \
    make && \
    make install && \
    cp -r * $BIN

## Install cmake
WORKDIR $SRC
RUN wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz && \
    tar xf cmake-3.2.2.tar.gz && \
    cd cmake-3.2.2 && \
    cp -r * ../. && \
    cd .. && \
    ./configure && \
    make && \
    make install && \
    cp -r * $BIN

RUN rm /usr/local/src/CMakeCache.txt
