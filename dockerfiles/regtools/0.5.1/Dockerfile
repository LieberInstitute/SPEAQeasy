FROM brcloudproject/ubuntu18.04-mpich

#  Not a binary, but place script where it's already on the path
COPY bed_to_juncs.py /usr/local/bin/bed_to_juncs.py

#  Install regtools
RUN wget https://github.com/griffithlab/regtools/archive/0.5.1.tar.gz -O regtools-0.5.1.tar.gz && \
    tar -xvf regtools-0.5.1.tar.gz && \
    cd regtools-0.5.1 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    mv regtools /usr/bin/
    
#  Install Python 2.7.9 (a bit more involved as package managers no longer
#  support Python < 3)
RUN wget https://www.python.org/ftp/python/2.7.9/Python-2.7.9.tgz && \
    tar xzvf Python-2.7.9.tgz && \
    chmod -R 755 Python-2.7.9 && \
    cd Python-2.7.9 && \
    ./configure && \
    make && \
    make install
