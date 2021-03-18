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

set -e

if [ "$1" == "docker" ]; then

    echo "User selected docker mode: nextflow will be installed, along with some test files."
    
    INSTALL_DIR=$(pwd)/Software
    mkdir -p $INSTALL_DIR
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
        -v $(pwd)/test:/test \
        $R_container \
        Rscript scripts/make_test_manifests.R -d $(pwd)
    
    #  Point to the original repository so that the "main" scripts can be
    #  trivially copied to share the pipeline
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_pipeline_sge.sh
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_pipeline_local.sh
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_pipeline_slurm.sh
        
elif [ "$1" == "local" ]; then

    echo -e "User selected local install: all software dependencies will be installed on this machine.\n\n"
    
    #  Verify java and python can be executed, since this is a pre-requisite
    if [ -x "$(command -v java)" ] && [ -x "$(command -v python3)" ]; then
        echo "Found Python 3 and a java runtime. Proceeding with the setup..."
  
        INSTALL_DIR=$(pwd)/Software
        mkdir -p $INSTALL_DIR
        cd $INSTALL_DIR
        
        #  Install nextflow (latest)
        wget -qO- https://get.nextflow.io | bash 
    
        #  bc (1.06.95)
        wget ftp://alpha.gnu.org/gnu/bc/bc-1.06.95.tar.bz2
        tar -xjf bc-1.06.95.tar.bz2
        cd bc-1.06.95
        ./configure --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        
        #  bcftools (1.10.2)  -------------------------------------------------------------
    
        wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 -O bcftools.tar.bz2
        tar -xjf bcftools.tar.bz2
        cd bcftools-1.10.2
        ./configure prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
            
        #  cmake, which regtools and salmon need to build
        git clone https://gitlab.kitware.com/cmake/cmake.git
        cd cmake
        ./bootstrap --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        
        #  fastqc (0.11.8)  -------------------------------------------------------------
    
        wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
        unzip fastqc_v0.11.8.zip
        chmod -R 775 FastQC
      
        #  hisat2 (2.2.1)  -------------------------------------------------------------
    
        wget https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.tar.gz
        tar -xzf v2.2.1.tar.gz
        cd hisat2-2.2.1
        make
        cd $INSTALL_DIR
        
        #  htslib (1.10.2)  -------------------------------------------------------------
    
        wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 -O htslib.tar.bz2
        tar -xjf htslib.tar.bz2
        cd htslib-1.10.2
        ./configure prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        mv htslib-1.10.2 htslib # wiggletools expects a "plain" name from htslib
            
        #  HDF5 (1.10.5), needed for build of kallisto -------------------------------
        
        wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
        tar -xzf hdf5-1.10.5.tar.gz
        cd hdf5-1.10.5
        ./configure --disable-parallel --without-szlib --without-pthread --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
    
        #  kallisto (0.46.1)  -------------------------------------------------------------
                
        wget https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
        tar -xzf v0.46.1.tar.gz
        cd kallisto-0.46.1
        mkdir build
        cd build
        $INSTALL_DIR/bin/cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
        make
        make prefix=$INSTALL_DIR install
        cd $INSTALL_DIR
    
        #  R (3.6.1) ---------------------------------------------------------------------
        
        #  Install R
        wget http://cran.rstudio.com/src/base/R-3/R-3.6.1.tar.gz
        tar -xf R-3.6.1.tar.gz
        cd R-3.6.1
        ./configure --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
      
        #  Install packages that will be used by the pipeline
        ./R-3.6.1/bin/Rscript ../scripts/check_R_packages.R
    
        #  Create the test samples.manifest files
        cd ..
        Software/R-3.6.1/bin/Rscript scripts/make_test_manifests.R -d $(pwd)
        cd $INSTALL_DIR
    
        #  regtools (0.5.1)  -------------------------------------------------------------
    
        wget https://github.com/griffithlab/regtools/archive/0.5.1.tar.gz -O regtools-0.5.1.tar.gz
        tar -xf regtools-0.5.1.tar.gz
        cd regtools-0.5.1
        mkdir build
        cd build
        $INSTALL_DIR/bin/cmake ..
        make
        cd $INSTALL_DIR
        
        #  rseqc (3.0.1)  -------------------------------------------------------------

        #  Install for the user, because we are not guaranteed system-wide
        #  installation privileges
        python3 -m pip install --user RSeQC==3.0.1
        
        #  salmon (1.2.1)  -------------------------------------------------------------
            
        wget https://github.com/COMBINE-lab/salmon/archive/v1.2.1.tar.gz
        tar -xzf v1.2.1.tar.gz
        mkdir build
        cd build
        $INSTALL_DIR/bin/cmake -DFETCH_BOOST=TRUE ..
        make
        make install
        cd $INSTALL_DIR
      
        #  samtools (1.10)  -------------------------------------------------------------
    
        wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 -O samtools.tar.bz2
        tar -xjf samtools.tar.bz2
        cd samtools-1.10
        ./configure prefix=$INSTALL_DIR
        make && \
        make install
        cd $INSTALL_DIR
        
        #  subread/ featureCounts (2.0.0)  -------------------------------------------------------------
            
        #  This works on Linux systems but needs testing on MacOS, FreeBSD
        if [ $(uname -s) == "Linux" ]; then
            makefile_name="Makefile.Linux"
        elif [ $(uname -s) == "Darwin" ]; then
            makefile_name="Makefile.MacOS"
        else
            makefile_name="Makefile.FreeBSD"
        fi
        
        wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-source.tar.gz/download
        tar -zxf download
        cd subread-2.0.0-source/src
        make -f $makefile_name
        cd $INSTALL_DIR
          
        #  trimmomatic (0.39)  -------------------------------------------------------------
        
        wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
        unzip Trimmomatic-0.39.zip
        chmod -R 775 Trimmomatic-0.39
          
        #  wiggletools (1.2.1)  -------------------------------------------------------------
        
        ## Install gsl
        wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
        tar -xf gsl-latest.tar.gz
        cd gsl-2.6
        ./configure --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        mv gsl-2.6 gsl # wiggletools expects a "plain" name from gsl
            
        ## Install libBigWig
        git clone git@github.com:dpryan79/libBigWig.git
        cd libBigWig
        make prefix=$INSTALL_DIR install
        cd $INSTALL_DIR
        
        # wiggletools itself (note the modified Makefile)
        git clone git@github.com:Ensembl/WiggleTools.git
        ./R-3.6.1/bin/Rscript ../scripts/fix_makefile.R
        cd WiggleTools
        make prefix=$INSTALL_DIR
        cd $INSTALL_DIR
    
        #  wigToBigWig -----------------------------------------------------------
        
        if [ $(uname -s) == "Linux" ]; then
            wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
            mv wigToBigWig bin/
        elif [ $(uname -s) == "Darwin" ]; then
            wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/wigToBigWig
            mv wigToBigWig bin/
        else
            #  Attempt to build from source without assuming the OS
            wget http://hgdownload.cse.ucsc.edu/admin/exe/userApps.src.tgz
            tar -xzf userApps.src.tgz
            cd userApps
            make
            mv bin/wigToBigWig ../bin/wigToBigWig
            cd $INSTALL_DIR
        fi
        
          
        #  Clean up compressed files
        rm $INSTALL_DIR/*.tar.gz
        rm -f $INSTALL_DIR/*.tgz
        rm $INSTALL_DIR/*.bz2
        rm $INSTALL_DIR/*.zip
        rm download
        
        #  Fix any strict permissions which would not allow sharing software
        #  (and therefore the pipeline as a whole) with those in a group
        chmod 775 -R $INSTALL_DIR
        
        #  Point to the original repository so that the "main" scripts can be
        #  trivially copied to share the pipeline
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(dirname $(pwd))|" ../run_pipeline_sge.sh
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(dirname $(pwd))|" ../run_pipeline_local.sh
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(dirname $(pwd))|" ../run_pipeline_slurm.sh
    
    else #  Java or Python could not be found on the system
    
        #  Java?
        if ! [ -x "$(command -v java)" ]; then
            echo "A java runtime could not be found or accessed. Is it installed and on the PATH? You can install it by running 'apt install default-jre', which requires sudo/ root privileges."
        fi
    
        #  Python?
        if ! [ -x "$(command -v python3)" ]; then
            echo "Python 3 could not be found or executed. Please install it (this requires root privileges)- ask an admin if needed."
        fi
        
        echo -e "\nAfter installing the required software, rerun this script to finish the installation procedure."
    fi
    
else # neither "docker" nor "local" was chosen
    
    echo 'Error: please specify "docker" or "local" and rerun this script.'
    echo '    eg. bash install_software.sh "local"'
    exit 1
    
fi
