FROM libddocker/ubuntu16.04_base:latest

WORKDIR /usr/local/src

#  Download trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    chmod -R 755 Trimmomatic-0.39

# Add the executable and adapters to the path
ENV PATH=$PATH:/usr/local/src/Trimmomatic-0.39:/usr/local/src/Trimmomatic-0.39/adapters
