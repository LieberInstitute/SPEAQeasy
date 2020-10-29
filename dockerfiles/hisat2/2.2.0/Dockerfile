FROM brcloudproject/ubuntu18.04-mpich

#  Install some tools/libraries
RUN apt-get update
RUN apt-get install -y wget unzip gcc g++ make zlib1g-dev libssl-dev libffi-dev liblzma-dev libbz2-dev liblzo2-dev zip
    
#  Install python 3
RUN apt-get install -y python3 2to3 python3-pip

RUN cd /usr/bin && \
    ln -s python3 python
    
RUN wget https://github.com/DaehwanKimLab/hisat2/archive/v2.2.0.tar.gz && \
    tar -xzf v2.2.0.tar.gz && \
    cd hisat2-2.2.0 && \
    make && \
    2to3 -w *.py && \
    2to3 -w */*.py && \
    sed -i "s/gzip.open(\(.*\), 'r')/gzip.open(\1, 'rt')/" hisat2_read_statistics.py && \
    cp hisat2 /usr/local/bin/ && \
    cp hisat2-align* /usr/local/bin/ && \
    cp hisat2-build* /usr/local/bin/ && \
    cp hisat2-inspect* /usr/local/bin/ && \
    cp *.py /usr/local/bin/
