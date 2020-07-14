FROM ubuntu:16.04

#  Install some tools/libraries
RUN apt-get update
RUN apt-get install -y wget bzip2 gcc make zlib1g-dev libssl-dev libffi-dev liblzma-dev libncurses5-dev libbz2-dev libcurl3-dev

#  Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjf samtools.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install
