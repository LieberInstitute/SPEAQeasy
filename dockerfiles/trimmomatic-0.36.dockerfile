FROM libddocker/ubuntu16.04_base:latest
 
MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"
 
########### Update Core ###########
RUN apt-get update
 
## define the folders for installation
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

## Install Trimmomatic
WORKDIR $SRC
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \ 
    unzip Trimmomatic-0.36.zip && \
    chmod -R 755 /usr/local/src/Trimmomatic-0.36 && \
    cp /usr/local/src/Trimmomatic-0.36/trimmomatic-0.36.jar $BIN && \
    cp /usr/local/src/Trimmomatic-0.36/adapters/* $BIN
