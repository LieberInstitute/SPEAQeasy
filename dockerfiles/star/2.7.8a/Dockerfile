FROM libddocker/ubuntu16.04_base:latest

RUN wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz && \
        tar -xzf 2.7.8a.tar.gz && \
        cd STAR-2.7.8a/source && \
        make STAR && \
        cp ../bin/Linux_x86_64/STAR /usr/local/bin/
