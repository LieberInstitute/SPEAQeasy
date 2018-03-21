FROM libddocker/ubuntu16.04_base:latest
 
MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"
 
########### Update Core ###########
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y cmake-curses-gui

## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin
ENV LIB /usr/local/lib


## Install cmake
WORKDIR $SRC
RUN wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz && \
    tar xf cmake-3.2.2.tar.gz && \
    cd cmake-3.2.2 && \
    cp -r * ../. && \
    cd .. && \
    ./configure && \
    make && \
    make install && \
    cp -r * $BIN

RUN rm /usr/local/src/CMakeCache.txt

## Install Salmon
WORKDIR $SRC
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz && \ 
    tar -xvzf Salmon-0.9.1_linux_x86_64.tar.gz && \
    cd Salmon-latest_linux_x86_64 && \
    mkdir build && \
    cmake .. && \
    make && \
    make install && \
    cp bin/* $BIN && \
    cp lib/* $LIB
