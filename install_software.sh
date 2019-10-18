#!/bin/bash
#$ -cwd

#  IMPORTANT: run this script inside the repository directory to ensure software
#              locations are properly linked. See conf/command_path_long.config if
#              you are interested in manually configuring different software paths.

#  Replace this line with 'local_install=true' if you wish to run this pipeline locally (rather
#  than through a cluster/ resource manager).
local_install=false

#  --------------------------------------------------------
#  Users should not need to alter code below this point
#  --------------------------------------------------------

if [ -n `which java` ]; then
  echo "Found a java runtime. Proceeding with the setup..."
  
  #  Install nextflow (latest)
  wget -qO- https://get.nextflow.io | bash
  
  if [ "$local_install" = true ]; then

    #  Downloads and installs all pipeline dependencies. The pipeline can then be run locally
    #  without further configuration.
    
    echo -e "User selected local install: all software dependencies will be installed on this machine.\n\n"
    INSTALL_DIR=$(pwd)/Software
    mkdir $INSTALL_DIR
    cd $INSTALL_DIR
    
    #  To do: consider how to handle java, wiggletools
    
    #  bcftools (1.9)  -------------------------------------------------------------
    
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 -O bcftools.tar.bz2 && \
        tar -xjvf bcftools.tar.bz2 && \
        cd bcftools-1.9 && \
        ./configure prefix=$INSTALL_DIR && \
        make && \
        make install
        cd $INSTALL_DIR
        
    #  fastqc (0.11.5)  -------------------------------------------------------------
    
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
      unzip fastqc_v0.11.5.zip && \
      chmod -R 755 FastQC
      
    #  hisat2 (2.1.0)  -------------------------------------------------------------
    
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
      unzip hisat2-2.1.0-Linux_x86_64.zip && \
      chmod -R 755 hisat2-2.1.0
      
    #  htslib (1.9)  -------------------------------------------------------------
    
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 -O htslib-1.9.tar.bz2 && \
        tar -xjvf htslib-1.9.tar.bz2 && \
        cd htslib-1.9 && \
        ./configure prefix=$INSTALL_DIR && \
        make && \
        make install
        cd $INSTALL_DIR
        mv htslib-1.9 htslib # wiggletools expects a "plain" name from htslib
    
    #  kallisto (0.43.0)  -------------------------------------------------------------
    
    wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz && \
      tar xzvf kallisto_linux-v0.43.0.tar.gz && \
      chmod -R 755 kallisto_linux-v0.43.0
    
    #  R (3.6.1) ---------------------------------------------------------------------
    
    #  Install R
    wget http://cran.rstudio.com/src/base/R-3/R-3.6.1.tar.gz && \
      tar xvf R-3.6.1.tar.gz && \
      cd R-3.6.1 && \
      ./configure --prefix=$INSTALL_DIR && \
      make && \
      make install
      cd $INSTALL_DIR
      
    #  Install packages that will be used by the pipeline
    ./R-3.6.1/bin/Rscript ../scripts/check_R_packages.R
    
    #  regtools (0.3.0)  -------------------------------------------------------------
    
    #  cmake, which regtools needs to build
    wget https://github.com/Kitware/CMake/releases/download/v3.14.6/cmake-3.14.6-Linux-x86_64.tar.gz && \
      tar xzvf cmake-3.14.6-Linux-x86_64.tar.gz
    
    #  regtools itself
    wget https://github.com/griffithlab/regtools/archive/0.3.0.tar.gz -O regtools-0.3.0.tar.gz && \
        tar -xvf regtools-0.3.0.tar.gz &&\
        cd regtools-0.3.0 && \
        mkdir build && \
        cd build && \
        ../../cmake-3.14.6-Linux-x86_64/bin/cmake .. && \
        make
        cd $INSTALL_DIR
        
    #  rseqc (2.6.4)  -------------------------------------------------------------
    
    #  Python (2.7.9), which RSeQC uses
    wget https://www.python.org/ftp/python/2.7.9/Python-2.7.9.tgz && \
      tar xzvf Python-2.7.9.tgz && \
      chmod -R 755 Python-2.7.9 && \
      cd Python-2.7.9 && \
      ./configure prefix=$INSTALL_DIR && \
      make && \
      make install
      cd $INSTALL_DIR
      
    #  RSeQC itself, and manual configuration of file locations so dependencies can be found 
    $INSTALL_DIR/Python-2.7.9/python -m pip install --root=temp -I RSeQC==2.6.4
    mv $INSTALL_DIR/temp/$INSTALL_DIR/bin/* $INSTALL_DIR/bin/
    mv $INSTALL_DIR/temp/$INSTALL_DIR/lib/python2.7/site-packages/* $INSTALL_DIR/bin/
    rm -rf $INSTALL_DIR/temp
        
    #  salmon (0.8.2)  -------------------------------------------------------------
    
    wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz && \
      tar xzvf Salmon-0.8.2_linux_x86_64.tar.gz && \
      chmod -R 755 Salmon-0.8.2_linux_x86_64
      
    #  samtools (1.9)  -------------------------------------------------------------
    
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -O samtools.tar.bz2 && \
        tar -xjvf samtools.tar.bz2 && \
        cd samtools-1.9 && \
        ./configure prefix=$INSTALL_DIR && \
        make && \
        make install
        cd $INSTALL_DIR
        
    #  subread/ featureCounts (1.5.0-p3)  -------------------------------------------------------------
    
    wget https://sourceforge.net/projects/subread/files/subread-1.5.0-p3/subread-1.5.0-p3-Linux-x86_64.tar.gz/download && \
      tar xzvf download && \
      chmod -R 755 subread-1.5.0-p3-Linux-x86_64
      
    #  trimmomatic (0.36)  -------------------------------------------------------------
    
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
      unzip Trimmomatic-0.36.zip && \
      chmod -R 755 Trimmomatic-0.36
      
    #  wiggletools (1.2.1)  -------------------------------------------------------------
    
    ## Install gsl
    wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && \
        tar xvf gsl-latest.tar.gz && \
        cd gsl-2.6 && \
        ./configure --prefix=$INSTALL_DIR && \
        make && \
        make install
        cd $INSTALL_DIR
        
    ## Install libBigWig
    git clone git@github.com:dpryan79/libBigWig.git && \
        cd libBigWig && \
        make prefix=$INSTALL_DIR install
        cd $INSTALL_DIR
    
    # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_DIR/lib
    git clone git@github.com:Ensembl/WiggleTools.git && \
        ./R-3.6.1/bin/Rscript ../scripts/fix_makefile.R && \
        cd WiggleTools && \
        make
        ## Make the one of the edits Mark Miller described at https://lists.johnshopkins.edu/sympa/arc/bithelp/2019-09/msg00132.html
    
    #  wigToBigWig -----------------------------------------------------------
    
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig && \
      chmod 755 wigToBigWig
      
    #  Clean up compressed files
    rm $INSTALL_DIR/*.tar.gz
    rm $INSTALL_DIR/*.tgz
    rm $INSTALL_DIR/*.bz2
    rm $INSTALL_DIR/*.zip
    
  else
    #  Already installed nextflow; simply need to install R packages used in the pipeline.
    #  Here it is assumed 'Rscript' is on the PATH (R is installed and accessible).
    echo -e "User selected typical install (bioinformatics tools will not be installed locally).\n\n"
    
    Rscript ../scripts/check_R_packages.R
  fi
else
  #  Java could not be found on the system
  echo "A java runtime could not be found or accessed. Is it installed and on the PATH? You can install it by running 'apt install default-jre', which may require root/sudo privileges. After proper installation, rerun this script to properly set up this pipeline."
fi
