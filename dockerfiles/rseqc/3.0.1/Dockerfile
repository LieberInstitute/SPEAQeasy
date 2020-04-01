FROM ubuntu:16.04

#  Install some tools/libraries
RUN apt-get update
RUN apt-get install -y wget gcc make zlib1g-dev libssl-dev libffi-dev liblzma-dev

#  Install python 3.7.7
RUN wget https://www.python.org/ftp/python/3.7.7/Python-3.7.7.tgz && \
    tar -xzf Python-3.7.7.tgz && \
    cd Python-3.7.7 && \
    ./configure && \
    make && \
    make install && \
    cp python /usr/bin/

#  Install RSeQC 3.0.1
RUN python -m pip install --upgrade pip
RUN python -m pip install rseqc==3.0.1
