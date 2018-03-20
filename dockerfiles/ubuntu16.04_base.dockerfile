FROM ubuntu:16.04

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
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
    liblzma-dev
RUN apt-get clean
