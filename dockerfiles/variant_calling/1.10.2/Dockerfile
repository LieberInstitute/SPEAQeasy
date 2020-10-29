FROM libddocker/samtools:1.10

#  Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjf bcftools.tar.bz2 && \
    cd bcftools-1.10.2 && \
    ./configure && \
    make && \
    make install

#  Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjf htslib.tar.bz2 && \
    cd htslib-1.10.2 && \
    ./configure && \
    make && \
    make install
