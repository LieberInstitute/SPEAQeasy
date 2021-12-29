FROM brcloudproject/ubuntu18.04-mpich

#  Install regtools (gpertea fork: v0.5.33g)
RUN wget https://github.com/gpertea/regtools/archive/refs/tags/0.5.33g.tar.gz && \
    tar -xf 0.5.33g.tar.gz && \
    cd regtools-0.5.33g && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    mv regtools /usr/bin/
