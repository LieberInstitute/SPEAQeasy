FROM ubuntu:18.04

#  Install some tools/libraries
RUN apt-get update
RUN apt-get install -y wget unzip gcc g++ make zlib1g-dev libssl-dev libffi-dev liblzma-dev libbz2-dev liblzo2-dev zip
    
#  Install python 3
RUN apt-get install -y python3

RUN cd /usr/bin && \
    ln -s python3 python
    
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.tar.gz && \
    tar -xzf v2.2.1.tar.gz && \
    cd hisat2-2.2.1 && \
    make && \
    cp hisat2 /usr/local/bin/ && \
    cp hisat2-align* /usr/local/bin/ && \
    cp hisat2-build* /usr/local/bin/ && \
    cp hisat2-inspect* /usr/local/bin/ && \
    cp *.py /usr/local/bin/
