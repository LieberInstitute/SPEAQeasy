FROM libddocker/ubuntu16.04_base:latest

RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar xzvf kallisto_linux-v0.46.1.tar.gz && \
    chmod -R 755 kallisto && \
    cp kallisto/kallisto /usr/local/bin/
