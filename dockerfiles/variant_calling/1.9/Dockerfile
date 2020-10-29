FROM libddocker/samtools:1.9

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.9 && \
    make && \
    make prefix=/usr/local/bin install && \
    cp -r * /usr/local/bin/
