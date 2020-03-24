#!/bin/bash

#  IMPORTANT: run this script inside the repository directory to ensure software
#              locations are properly linked. See conf/command_path_long.config if
#              you are interested in manually configuring different software paths.
#
#  Usage:  bash install_software.sh [installation_type]
#
#    installation_type may be "docker" or "local":
#        docker:    user plans to run pipeline with docker to manage software dependencies
#        local:     user wishes to install all dependencies locally (regardless of whether
#                   the pipeline will be run on a cluster or with local resources) 


#  This is the docker image to be used for execution of R via docker (with docker mode)
R_container="libddocker/r_3.6.1_bioc"

#  --------------------------------------------------------
#  Users should not need to alter code below this point
#  --------------------------------------------------------

if [ "$1" == "docker" ]; then

    echo "User selected docker mode: nextflow will be installed, along with some test files."
    
    INSTALL_DIR=$(pwd)/Software
    mkdir $INSTALL_DIR
    cd $INSTALL_DIR
    
    #  Install nextflow (latest)
    wget -qO- https://get.nextflow.io | bash
    cd ..
        
    ###########################################################################
    #  Create the samples.manifest files for test directories
    ###########################################################################
    
    #  Grab the container
    docker pull $R_container
    
    docker run \
        -v $(pwd)/scripts:/scripts/ \
        $R_container \
        Rscript scripts/make_test_manifests.R
        
elif [ "$1" == "local" ]; then

    echo -e "User selected local install: all software dependencies will be installed on this machine.\n\n"
    
    #  Verify java and python can be executed, since this is a pre-requisite
    if [ -x "$(command -v java)" ] && [ -x "$(command -v python)" ]; then
        #  Ensure python version is python 3
        if [ ! "$(python -V | cut -d " " -f 2 | cut -d "." -f 1)" == "3" ]; then
            echo "Python is installed, but not python 3 in particular. Python 3 is a prerequisite for SPEQeasy."
            exit 1
        fi

        echo "Found Python 3 and a java runtime. Proceeding with the setup..."
  
        INSTALL_DIR=$(pwd)/Software
        mkdir $INSTALL_DIR
        cd $INSTALL_DIR
  
        #  Install nextflow (latest)
        wget -qO- https://get.nextflow.io | bash 
    
        #  bc (1.06.95)
        wget ftp://alpha.gnu.org/gnu/bc/bc-1.06.95.tar.bz2 && \
            tar -xjvf bc-1.06.95.tar.bz2 && \
            cd bc-1.06.95 && \
            ./configure --prefix=$INSTALL_DIR
            make
            make install
        
        #  bcftools (1.9)  -------------------------------------------------------------
    
        wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 -O bcftools.tar.bz2 && \
            tar -xjvf bcftools.tar.bz2 && \
            cd bcftools-1.9 && \
            ./configure prefix=$INSTALL_DIR && \
            make && \
            make install
            cd $INSTALL_DIR
        
        #  fastqc (0.11.8)  -------------------------------------------------------------
    
        wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip && \
            unzip fastqc_v0.11.8.zip && \
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
    
        #  kallisto (0.46.1)  -------------------------------------------------------------
    
        wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
            tar xzvf kallisto_linux-v0.46.1.tar.gz && \
            chmod -R 755 kallisto
    
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
    
        #  Create the test samples.manifest files
        ./R-3.6.1/bin/Rscript ../scripts/make_test_manifests.R
    
        #  regtools (0.5.1)  -------------------------------------------------------------
    
        #  cmake, which regtools needs to build
        wget https://github.com/Kitware/CMake/releases/download/v3.14.6/cmake-3.14.6-Linux-x86_64.tar.gz && \
            tar xzvf cmake-3.14.6-Linux-x86_64.tar.gz
    
        #  regtools itself
        wget https://github.com/griffithlab/regtools/archive/0.5.1.tar.gz -O regtools-0.5.1.tar.gz && \
            tar -xvf regtools-0.5.1.tar.gz &&\
            cd regtools-0.5.1 && \
            mkdir build && \
            cd build && \
            ../../cmake-3.14.6-Linux-x86_64/bin/cmake .. && \
            make
            cd $INSTALL_DIR
        
        #  rseqc (3.0.1)  -------------------------------------------------------------

        python -m pip install --user RSeQC==3.0.1
        
        #  salmon (1.0.0)  -------------------------------------------------------------
    
        wget https://github.com/COMBINE-lab/salmon/releases/download/v1.0.0/Salmon-1.0.0_linux_x86_64.tar.gz && \
            tar xzvf Salmon-1.0.0_*.tar.gz && \
            rm Salmon-1.0.0_*.tar.gz && \
            chmod -R 755 salmon-*
      
        #  samtools (1.9)  -------------------------------------------------------------
    
        wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -O samtools.tar.bz2 && \
            tar -xjvf samtools.tar.bz2 && \
            cd samtools-1.9 && \
            ./configure prefix=$INSTALL_DIR && \
            make && \
            make install
            cd $INSTALL_DIR
        
        #  subread/ featureCounts (2.0.0)  -------------------------------------------------------------
        
        wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz/download && \
            tar xzvf download && \
            chmod -R 755 subread-2.0.0-Linux-x86_64
          
        #  trimmomatic (0.39)  -------------------------------------------------------------
        
        wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
            unzip Trimmomatic-0.39.zip && \
            chmod -R 755 Trimmomatic-0.39
          
        #  wiggletools (1.2.1)  -------------------------------------------------------------
        
        ## Install gsl
        wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz && \
            tar xvf gsl-latest.tar.gz && \
            cd gsl-2.6 && \
            ./configure --prefix=$INSTALL_DIR && \
            make && \
            make install
            cd $INSTALL_DIR
            mv gsl-2.6 gsl # wiggletools expects a "plain" name from gsl
            
        ## Install libBigWig
        git clone git@github.com:dpryan79/libBigWig.git && \
            cd libBigWig && \
            make prefix=$INSTALL_DIR install
            cd $INSTALL_DIR
        
        # wiggletools itself (note the modified Makefile)
        git clone git@github.com:Ensembl/WiggleTools.git && \
            ./R-3.6.1/bin/Rscript ../scripts/fix_makefile.R && \
            cd WiggleTools && \
            make prefix=$INSTALL_DIR
            cd $INSTALL_DIR
    
        #  wigToBigWig -----------------------------------------------------------
        
        wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig && \
            chmod 755 wigToBigWig
          
        #  Clean up compressed files
        rm $INSTALL_DIR/*.tar.gz
        rm $INSTALL_DIR/*.bz2
        rm $INSTALL_DIR/*.zip
        rm download
    
    else #  Java or Python could not be found on the system
    
        #  Java?
        if ! [ -x "$(command -v java)" ]; then
            echo "A java runtime could not be found or accessed. Is it installed and on the PATH? You can install it by running 'apt install default-jre', which requires sudo/ root privileges."
        fi
    
        #  Python?
        if ! [ -x "$(command -v python)" ]; then
            echo "Python 3 could not be found or executed. Please install it (this requires root privileges)- ask an admin if needed."
        fi
    fi
  
    echo -e "\nAfter installing the required software, rerun this script to finish the installation procedure."
    
else # neither "docker" nor "local" was chosen
    
    echo 'Error: please specify "docker" or "local" and rerun this script.'
    echo '    eg. bash install_software.sh "local"'
    exit 1
    
fi
