FROM libddocker/ubuntu16.04_base:latest

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.14.1/Salmon-0.14.1_linux_x86_64.tar.gz && \
    tar xzvf Salmon-0.14.1_*.tar.gz && \
    rm Salmon-0.14.1_*.tar.gz && \
    chmod -R 755 salmon-* && \
    mv salmon-*/bin/salmon /usr/local/bin/ && \
    mv salmon-*/lib/* /usr/local/lib/
