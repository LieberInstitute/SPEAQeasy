FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip && \
    chmod -R 755 hisat2-2.1.0 && \
    cp hisat2-2.1.0/hisat2 /usr/local/bin/ && \
    cp hisat2-2.1.0/hisat2-align* /usr/local/bin/ && \
    cp hisat2-2.1.0/hisat2-build* /usr/local/bin/ && \
    cp hisat2-2.1.0/hisat2-inspect* /usr/local/bin/
