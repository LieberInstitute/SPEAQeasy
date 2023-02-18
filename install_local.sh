#!/bin/bash

#  Usage:  bash install_local.sh
#
#   This script only installs auxiliary programs that are not already found 
#   in PATH
#   Recent versions are assumed - note that regtools MUST be the one from
#   https://github.com/gpertea/regtools

set -e

REPO_NAME=$(basename $(pwd))
#    echo "Making sure the repository is clean and ready for installation..."
#    
#    #rm -r Software
#    rm -f test/*/*/*/samples.manifest
#    git checkout run_pipeline_*.sh
#    git checkout nextflow.config
#    git checkout conf/*.config
#    

error_message() {
    if [[ $(uname -s) == "Darwin" ]] && [[ "$1" == "local" ]]; then
        echo -e "\nNote that 'local' installation is not officially supported on Mac OS! Please use the 'docker' mode instead if possible."
    fi
    
    echo -e "\nPlease check http://research.libd.org/SPEAQeasy/setup-details.html#troubleshooting for information about resolving installation-related issues. Note that the 'Software' directory should be deleted before re-attempting installation."
}

trap "error_message" ERR

    #  Verify java and python can be executed, since this is a pre-requisite
    if [[ -x "$(command -v java)" && -x "$(command -v python3)" && -x "$(command -v R)" ]]; then
        echo "Found Python 3, java and R. Proceeding with the setup..."
  
        INSTALL_DIR=$(pwd)/Software
        mkdir -p $INSTALL_DIR
        cd $INSTALL_DIR
        
        #  Install nextflow (20.01.0)
        #curl -L https://github.com/nextflow-io/nextflow/releases/download/v20.01.0/nextflow | bash
    
        #  bc (1.06.95)
        #curl -O ftp://alpha.gnu.org/gnu/bc/bc-1.06.95.tar.bz2
        #tar -xjf bc-1.06.95.tar.bz2
        #cd bc-1.06.95
        #./configure --prefix=$INSTALL_DIR
        #make
        #make install
        #cd $INSTALL_DIR
        
        #  bcftools (1.10.2)  -------------------------------------------------------------
    
        #curl -o bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
        #tar -xjf bcftools.tar.bz2
        #cd bcftools-1.10.2
        #./configure prefix=$INSTALL_DIR
        #make
        #make install
        #cd $INSTALL_DIR
            
        #  cmake, which regtools and salmon need to build
        #git clone https://gitlab.kitware.com/cmake/cmake.git
        #cd cmake
        #./bootstrap --prefix=$INSTALL_DIR
        #make
        #make install
        #cd $INSTALL_DIR
        
        #  fastqc (0.11.8)  -------------------------------------------------------------
        if ! command -v fastqc &> /dev/null ; then
          curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
          unzip fastqc_v0.11.8.zip
          chmod -R 775 FastQC
        fi
      
        #  hisat2 (2.2.1)  -------------------------------------------------------------
       if ! command -v hisat2 &> /dev/null ; then
        curl -O https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.tar.gz
        tar -xzf v2.2.1.tar.gz
        cd hisat2-2.2.1
        make
        cd $INSTALL_DIR
       fi
        #  htslib (1.10.2)  -------------------------------------------------------------
    
        #curl -o htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
        #tar -xjf htslib.tar.bz2
        #cd htslib-1.10.2
        #./configure prefix=$INSTALL_DIR
        #make
        #make install
        #cd $INSTALL_DIR
        #mv htslib-1.10.2 htslib # wiggletools expects a "plain" name from htslib
            
        #  HDF5 (1.10.5), needed for build of kallisto -------------------------------
        
        #curl -O https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
        #tar -xzf hdf5-1.10.5.tar.gz
        #cd hdf5-1.10.5
        #./configure --disable-parallel --without-szlib --without-pthread --prefix=$INSTALL_DIR
        #make
        #make install
        #cd $INSTALL_DIR
    
        #  kallisto (0.46.1)  -------------------------------------------------------------
                
        #curl -O https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
        #tar -xzf v0.46.1.tar.gz
        #cd kallisto-0.46.1
        #mkdir build
        #cd build
        #$INSTALL_DIR/bin/cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
        #make
        #make prefix=$INSTALL_DIR install
        #cd $INSTALL_DIR
    
        #  R (4.1.2) ---------------------------------------------------------------------
        
        #  Install R
        #curl -O http://cran.rstudio.com/src/base/R-4/R-4.1.2.tar.gz
        #tar -xzf R-4.1.2.tar.gz
        #cd R-4.1.2
        #./configure --prefix=$INSTALL_DIR --with-x=no
        #make
        #make install
        #cd $INSTALL_DIR
      
        #  Install packages that will be used by the pipeline
        Rscript ../scripts/check_R_packages_local.R
        
        BASE_DIR=$(dirname $INSTALL_DIR)
        cd $BASE_DIR
        
        #  Signal to load ordinary R packages with 'checkpoint' in each R script
        #sed -i "1i #  Added during installation\nlibrary('checkpoint')\ncheckpoint('2021-10-01',\n    project_dir = '$BASE_DIR/scripts/r_packages',\n    checkpoint_location = '$BASE_DIR/Software'\n)\n" scripts/*.R
    
        #  Create the test samples.manifest files
        #Rscript scripts/make_test_manifests.R -d $(pwd)
        cd $INSTALL_DIR
    
        #  regtools (gpertea fork: 0.5.33g)  ----------------------------------
        if ! command -v regtools &> /dev/null ; then
          curl -O https://github.com/gpertea/regtools/archive/refs/tags/0.5.33g.tar.gz
          tar -xf 0.5.33g.tar.gz
          cd regtools-0.5.33g
          mkdir build
          cd build
          $INSTALL_DIR/bin/cmake ..
          make
          cd $INSTALL_DIR
        fi
        #  rseqc (3.0.1)  -------------------------------------------------------------

        #  Install for the user, because we are not guaranteed system-wide
        #  installation privileges
        python3 -m pip install --user RSeQC==3.0.1
        
        #  salmon (1.2.1)  -------------------------------------------------------------
            
        #if ! command -v salmon &> /dev/null ; then
        #  #wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz
        #  echo "Please install salmon!"
        #  exit 1
        #fi
        #  samtools (1.10)  -------------------------------------------------------------
    
#         curl -o samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
#         tar -xjf samtools.tar.bz2
#         cd samtools-1.10
#         ./configure prefix=$INSTALL_DIR
#         make && \
#         make install
#         cd $INSTALL_DIR
        
        #  STAR (2.7.8a)  ---------------------------------------------------------------
            
#         curl -O https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
#         tar -xzf 2.7.8a.tar.gz
#         cd STAR-2.7.8a/source
#         make STAR
#         cd $INSTALL_DIR
#         
#         #  Copy the binary to a fixed location, regardless of OS
#         if [ $(uname -s) == "Linux" ]; then
#             cp STAR-2.7.8a/bin/Linux_x86_64/STAR bin/
#         elif [ $(uname -s) == "Darwin" ]; then
#             cp STAR-2.7.8a/bin/MacOSX_x86_64/STAR bin/
#         else
#             echo "Non-docker installation of the STAR aligner is not supported on this operating system. Consider instead setting up SPEAQeasy for use with docker, by running:"
#             echo '    bash install_software.sh "docker"'
#             exit 1
#         fi
#         
#         #  subread/ featureCounts (2.0.0)  -------------------------------------------------------------
#             
#         #  This works on Linux systems but needs testing on MacOS, FreeBSD
#         if [ $(uname -s) == "Linux" ]; then
#             makefile_name="Makefile.Linux"
#         elif [ $(uname -s) == "Darwin" ]; then
#             makefile_name="Makefile.MacOS"
#         else
#             makefile_name="Makefile.FreeBSD"
#         fi
        
#         curl -O https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-source.tar.gz/download
#         tar -zxf download
#         cd subread-2.0.0-source/src
#         make -f $makefile_name
#         cd $INSTALL_DIR
          
        #  trimmomatic (0.39)  -------------------------------------------------------------
        
        curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
        unzip Trimmomatic-0.39.zip
        chmod -R 775 Trimmomatic-0.39
          
        #  wiggletools (1.2.1)  -------------------------------------------------------------
        
        ## Install gsl
#         curl -O https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
#         tar -xf gsl-2.6.tar.gz
#         cd gsl-2.6
#         ./configure --prefix=$INSTALL_DIR
#         make
#         make install
#         cd $INSTALL_DIR
#         mv gsl-2.6 gsl # wiggletools expects a "plain" name from gsl
#             
#         ## Install libBigWig
#         curl -O https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.6.tar.gz
#         tar -xzf 0.4.6.tar.gz
#         mv libBigWig-0.4.6 libBigWig
#         cd libBigWig
#         make prefix=$INSTALL_DIR install
#         cd $INSTALL_DIR
#         
#         # wiggletools itself (note the modified Makefile)
#         curl -O https://github.com/Ensembl/WiggleTools/archive/refs/tags/v1.2.1.tar.gz
#         tar -xzf v1.2.1.tar.gz
#         mv WiggleTools-1.2.1 WiggleTools
#         ./R-3.6.1/bin/Rscript ../scripts/fix_makefile.R
#         cd WiggleTools
#         make prefix=$INSTALL_DIR
#         cd $INSTALL_DIR
    
        #  wigToBigWig -----------------------------------------------------------
        
        #if [ $(uname -s) == "Linux" ]; then
        #    curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
        #    mv wigToBigWig bin/
        #elif [ $(uname -s) == "Darwin" ]; then
        #    curl -O http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/wigToBigWig
        #    mv wigToBigWig bin/
        #else
        #    #  Attempt to build from source without assuming the OS
        #    curl -O http://hgdownload.cse.ucsc.edu/admin/exe/userApps.src.tgz
        #    tar -xzf userApps.src.tgz
        #    cd userApps
        #    make
        #    mv bin/wigToBigWig ../bin/wigToBigWig
        #    cd $INSTALL_DIR
        #fi
        
          
        #  Clean up compressed files
        rm $INSTALL_DIR/*.tar.gz
        rm -f $INSTALL_DIR/*.tgz
        rm $INSTALL_DIR/*.bz2
        rm $INSTALL_DIR/*.zip
        rm -f download
        
        #  Fix any strict permissions which would not allow sharing software
        #  (and therefore the pipeline as a whole) with those in a group
        chmod 775 -R $INSTALL_DIR
        
        #  Point to the original repository so that the "main" scripts can be
        #  trivially copied to share the pipeline
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(dirname $(pwd))|" ../run_pipeline_sge.sh
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(dirname $(pwd))|" ../run_pipeline_local.sh
        sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(dirname $(pwd))|" ../run_pipeline_slurm.sh
    fi
