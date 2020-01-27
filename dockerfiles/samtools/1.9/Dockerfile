FROM libddocker/ubuntu16.04_base:latest

WORKDIR $SRC

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-1.9 && \
    make && \
    make prefix=/usr/local/bin install && \
    cp -r * $BIN
