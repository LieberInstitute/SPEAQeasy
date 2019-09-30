INSTALL_DIR=$(pwd)

#  To do: consider how to handle java, R 3.6, ucsctools/wigToBigWig

#  bcftools (1.9)

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.9 && \
    ./configure prefix=$INSTALL_DIR && \
    make && \
    make install
    
#  fastqc (0.11.5)

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
  unzip fastqc_v0.11.5.zip && \
  chmod -R 755 FastQC
  
#  hisat2 (2.1.0)

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
  unzip hisat2-2.1.0-Linux_x86_64.zip && \
  chmod -R 755 hisat2-2.1.0
  
#  htslib (1.9)

wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 -O htslib-1.9.tar.bz2 && \
    tar -xjvf htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure prefix=$INSTALL_DIR && \
    make && \
    make install

#  kallisto (0.43.0)

wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz && \
  tar xzvf kallisto_linux-v0.43.0.tar.gz && \
  chmod -R 755 kallisto_linux-v0.43.0
  
#  nextflow (latest)

wget -qO- https://get.nextflow.io | bash

#  python (2.7.9)

pip install python==2.7.9

#  regtools (0.3.0)

wget https://github.com/griffithlab/regtools/archive/0.3.0.tar.gz -O regtools-0.3.0.tar.gz && \
    tar -xvf regtools-0.3.0.tar.gz &&\
    cd regtools-0.3.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make
    
#  rseqc (2.6.4)

wget https://downloads.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz && \
  tar zxf RSeQC-2.6.4.tar.gz && \
  cd RSeQC-2.6.4 && \
  python setup.py install --root=bin
    
#  salmon (0.8.2)

wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz && \
  tar xzvf Salmon-0.8.2_linux_x86_64.tar.gz && \
  chmod -R 755 Salmon-0.8.2_linux_x86_64
  
#  samtools (1.9)

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure prefix=$INSTALL_DIR && \
    make && \
    make install
    
#  subread/ featureCounts (1.5.0-p3)

wget https://sourceforge.net/projects/subread/files/subread-1.5.0-p3/subread-1.5.0-p3-Linux-x86_64.tar.gz/download && \
  tar xzvf download && \
  chmod -R 755 subread-1.5.0-p3-Linux-x86_64
  
#  trimmomatic (0.36)

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
  unzip Trimmomatic-0.36.zip && \
  chmod -R 755 Trimmomatic-0.36
  
#  wiggletools (1.2.1)

## Install gsl
wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && \
    tar xvf gsl-latest.tar.gz && \
    cd gsl-2.6 && \
    ./configure --prefix=$INSTALL_DIR && \
    make && \
    make install
    
## Install libBigWig
git clone git@github.com:dpryan79/libBigWig.git && \
    cd libBigWig && \
    make prefix=$INSTALL_DIR install
    
git clone git@github.com:Ensembl/WiggleTools.git && \
    cd WiggleTools && \
    ## Make the one of the edits Mark Miller described at https://lists.johnshopkins.edu/sympa/arc/bithelp/2019-09/msg00132.html
    ## Change LIBS from:
    # LIBS= -lwiggletools -l:libBigWig.a -lcurl -l:libhts.a -lgsl  -lgslcblas -lz -lpthread -lm
    ## to:
    # LIBS= -lwiggletools -l:libBigWig.a -lcurl -l:libhts.a -lgsl  -lgslcblas -lz -lpthread -lm -lcrypto -llzma -lbz2
    make
    ## Copy to bin folder with other stuff
    cp bin/* ../bin/



