FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
    unzip fastqc_v0.11.8.zip && \
    chmod -R 755 FastQC && \
    ln -s /usr/local/src/FastQC/fastqc /usr/local/bin
