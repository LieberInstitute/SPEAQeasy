FROM libddocker/ubuntu16.04_base:latest

RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-source.tar.gz/download && \
    tar -zxvf download && \
    cd subread-2.0.0-source/src && \
    make -f Makefile.Linux && \
    cp -R ../bin/* /usr/local/bin/