FROM libddocker/ubuntu16.04_base:latest

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.2.1/Salmon-1.2.1_linux_x86_64.tar.gz && \
    tar xzvf Salmon-1.2.1_*.tar.gz && \
    rm Salmon-1.2.1_*.tar.gz && \
    chmod -R 775 salmon-* && \
    mv salmon-*/bin/salmon /usr/local/bin/ && \
    mv salmon-*/lib/* /usr/local/lib/
