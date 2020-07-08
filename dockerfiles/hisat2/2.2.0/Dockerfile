FROM ubuntu:16.04

#  Install some tools/libraries
RUN apt-get update
RUN apt-get install -y wget unzip gcc g++ make zlib1g-dev libssl-dev libffi-dev liblzma-dev

#  Install python 3.7.7
RUN wget https://www.python.org/ftp/python/3.7.7/Python-3.7.7.tgz && \
    tar -xzf Python-3.7.7.tgz && \
    cd Python-3.7.7 && \
    ./configure && \
    make && \
    make install && \
    cp python /usr/bin/ && \
    cd ..
    
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/v2.2.0.tar.gz && \
    tar -xzf v2.2.0.tar.gz && \
    cd hisat2-2.2.0 && \
    make && \
    cp hisat2 /usr/local/bin/ && \
    cp hisat2-align* /usr/local/bin/ && \
    cp hisat2-build* /usr/local/bin/ && \
    cp hisat2-inspect* /usr/local/bin/
