#!/bin/bash

#  Usage:  bash install_software.sh [installation_type]
#
#    installation_type may be "docker", "singularity", "local", or "jhpce":
#        docker:      user plans to run pipeline with docker to manage software
#                     dependencies
#        singularity: user plans to run pipeline with singularity to manage
#                     software
#        local:       user wishes to install all dependencies locally
#                     (regardless of whether the pipeline will be run on a
#                     cluster or with local resources) 
#        jhpce:       user is setting up SPEAQeasy on the JHPCE cluster


#  This is the docker image to be used for execution of R via docker or
#  singularity, if applicable
R_container="libddocker/bioc_kallisto:3.13"

#  --------------------------------------------------------
#  Users should not need to alter code below this point
#  --------------------------------------------------------

set -e

REPO_NAME="SPEAQeasy"

if [ "$(basename $(pwd))" == "$REPO_NAME" ]; then

    echo "Making sure the repository is clean and ready for installation..."
    
    #rm -r Software
    rm -f test/*/*/*/samples.manifest
    git checkout run_pipeline_*.sh
    git checkout nextflow.config
    git checkout conf/*.config
    
else

    echo "Please only invoke the script from directly inside the '$REPO_NAME' directory!"
    exit 1
    
fi


error_message() {
    if [[ $(uname -s) == "Darwin" ]] && [[ "$1" == "local" ]]; then
        echo -e "\nNote that 'local' installation is not officially supported on Mac OS! Please use the 'docker' mode instead if possible."
    fi
    
    echo -e "\nPlease check http://research.libd.org/SPEAQeasy/setup-details.html#troubleshooting for information about resolving installation-related issues. Note that the 'Software' directory should be deleted before re-attempting installation."
}

trap "error_message" ERR

if [[ "$1" == "docker" || "$1" == "singularity" ]]; then

    echo "User selected $1 mode: nextflow will be installed, along with some test files."
    
    INSTALL_DIR=$(pwd)/Software
    mkdir -p $INSTALL_DIR
    cd $INSTALL_DIR
    
    #  Install nextflow (latest)
    curl -s https://get.nextflow.io | bash
    cd ..
        
    ###########################################################################
    #  Create the samples.manifest files for test directories
    ###########################################################################
    
    if [[ "$1" == "docker" ]]; then
        if [[ "$2" == "sudo" ]]; then
            command="sudo docker run"
        else
            command="docker run"
        fi
        
        $command \
            -u $(id -u):$(id -g) \
            -v $(pwd)/scripts:/usr/local/src/scripts/ \
            -v $(pwd)/test:/usr/local/src/test \
            $R_container \
            Rscript /usr/local/src/scripts/make_test_manifests.R -d $(pwd)
    else # using singularity
        if [[ "$(singularity version | cut -d '.' -f 1)" -lt 3 ]]; then
            echo "Singularity >= 3 is required."
            exit 1
        fi
        
        echo "Pulling images in advance to have them cached for the first execution..."
        
        #  Pull images in advance, since it seems to use very large amounts of
        #  memory to build the '.sif' file from each docker image (we don't
        #  want to allocate large amounts of memory in each process just for
        #  this purpose)
        rm -rf dockerfiles/singularity_cache
        mkdir dockerfiles/singularity_cache
        images=$(grep 'container = ' conf/singularity.config | tr -d " |'" | cut -d '=' -f 2 | sort -u)
        for image in $images; do
            image_name=$(echo $image | sed 's/[:\/]/-/g').img
            singularity pull dockerfiles/singularity_cache/$image_name docker://$image
        done
        
        singularity exec \
            --pwd /usr/local/src \
            -B $(pwd)/scripts:/usr/local/src/scripts/ \
            -B $(pwd)/test:/usr/local/src/test \
            docker://$R_container \
            Rscript /usr/local/src/scripts/make_test_manifests.R -d $(pwd)
        
        #  Set modules used correctly for JHPCE users
        sed -i "/module = '.*\/.*'/d" conf/jhpce.config
        sed -i "s|cache = 'lenient'|cache = 'lenient'\n    module = 'singularity/3.6.0'|" conf/jhpce.config
        sed -i "s|module load nextflow|module load nextflow\nmodule load singularity/3.6.0|" run_pipeline_jhpce.sh
    fi
    
    #  Add docker/ singularity configuration to each config profile in
    #  'nextflow.config'; use short command paths instead of long
    sed -i -r "s:includeConfig 'conf/(local|sge|slurm|jhpce)\.config':includeConfig 'conf/\1.config'\n        includeConfig 'conf/$1.config':" nextflow.config
    sed -i "s|includeConfig 'conf/command_paths_long.config'|includeConfig 'conf/command_paths_short.config'|" nextflow.config
    
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
        curl -s https://get.nextflow.io | bash
    
        #  bc (1.06.95)
        curl -O ftp://alpha.gnu.org/gnu/bc/bc-1.06.95.tar.bz2
        tar -xjf bc-1.06.95.tar.bz2
        cd bc-1.06.95
        ./configure --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        
        #  bcftools (1.10.2)  -------------------------------------------------------------
    
        curl -o bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
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
    
        curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
        unzip fastqc_v0.11.8.zip
        chmod -R 775 FastQC
      
        #  hisat2 (2.2.1)  -------------------------------------------------------------
    
        curl -O https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.tar.gz
        tar -xzf v2.2.1.tar.gz
        cd hisat2-2.2.1
        make
        cd $INSTALL_DIR
        
        #  htslib (1.10.2)  -------------------------------------------------------------
    
        curl -o htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
        tar -xjf htslib.tar.bz2
        cd htslib-1.10.2
        ./configure prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        mv htslib-1.10.2 htslib # wiggletools expects a "plain" name from htslib
            
        #  HDF5 (1.10.5), needed for build of kallisto -------------------------------
        
        curl -O https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
        tar -xzf hdf5-1.10.5.tar.gz
        cd hdf5-1.10.5
        ./configure --disable-parallel --without-szlib --without-pthread --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
    
        #  kallisto (0.46.1)  -------------------------------------------------------------
                
        curl -O https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
        tar -xzf v0.46.1.tar.gz
        cd kallisto-0.46.1
        mkdir build
        cd build
        $INSTALL_DIR/bin/cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
        make
        make prefix=$INSTALL_DIR install
        cd $INSTALL_DIR
    
        #  R (4.1.0) ---------------------------------------------------------------------
        
        #  Install R
        curl -O http://cran.rstudio.com/src/base/R-4/R-4.1.0.tar.gz
        tar -xzf R-4.1.0.tar.gz
        cd R-4.1.0
        ./configure --prefix=$INSTALL_DIR --with-x=no
        make
        make install
        cd $INSTALL_DIR
      
        #  Install packages that will be used by the pipeline
        ./R-4.1.0/bin/Rscript ../scripts/check_R_packages.R
    
        #  Create the test samples.manifest files
        cd ..
        Software/R-4.1.0/bin/Rscript scripts/make_test_manifests.R -d $(pwd)
        cd $INSTALL_DIR
    
        #  regtools (0.5.1)  -------------------------------------------------------------
    
        curl -o regtools-0.5.1.tar.gz https://github.com/griffithlab/regtools/archive/0.5.1.tar.gz 
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
            
        curl -o salmon_v1.2.1.tar.gz https://github.com/COMBINE-lab/salmon/archive/v1.2.1.tar.gz
        tar -xzf salmon_v1.2.1.tar.gz
        cd salmon-1.2.1
        mkdir build
        cd build
        $INSTALL_DIR/bin/cmake -DFETCH_BOOST=TRUE ..
        make
        make install
        cd $INSTALL_DIR
      
        #  samtools (1.10)  -------------------------------------------------------------
    
        curl -o samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
        tar -xjf samtools.tar.bz2
        cd samtools-1.10
        ./configure prefix=$INSTALL_DIR
        make && \
        make install
        cd $INSTALL_DIR
        
        #  STAR (2.7.8a)  ---------------------------------------------------------------
            
        curl -O https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
        tar -xzf 2.7.8a.tar.gz
        cd STAR-2.7.8a/source
        make STAR
        cd $INSTALL_DIR
        
        #  Copy the binary to a fixed location, regardless of OS
        if [ $(uname -s) == "Linux" ]; then
            cp STAR-2.7.8a/bin/Linux_x86_64/STAR bin/
        elif [ $(uname -s) == "Darwin" ]; then
            cp STAR-2.7.8a/bin/MacOSX_x86_64/STAR bin/
        else
            echo "Non-docker installation of the STAR aligner is not supported on this operating system. Consider instead setting up SPEAQeasy for use with docker, by running:"
            echo '    bash install_software.sh "docker"'
            exit 1
        fi
        
        #  subread/ featureCounts (2.0.0)  -------------------------------------------------------------
            
        #  This works on Linux systems but needs testing on MacOS, FreeBSD
        if [ $(uname -s) == "Linux" ]; then
            makefile_name="Makefile.Linux"
        elif [ $(uname -s) == "Darwin" ]; then
            makefile_name="Makefile.MacOS"
        else
            makefile_name="Makefile.FreeBSD"
        fi
        
        curl -O https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-source.tar.gz/download
        tar -zxf download
        cd subread-2.0.0-source/src
        make -f $makefile_name
        cd $INSTALL_DIR
          
        #  trimmomatic (0.39)  -------------------------------------------------------------
        
        curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
        unzip Trimmomatic-0.39.zip
        chmod -R 775 Trimmomatic-0.39
          
        #  wiggletools (1.2.1)  -------------------------------------------------------------
        
        ## Install gsl
        curl -O https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
        tar -xf gsl-2.6.tar.gz
        cd gsl-2.6
        ./configure --prefix=$INSTALL_DIR
        make
        make install
        cd $INSTALL_DIR
        mv gsl-2.6 gsl # wiggletools expects a "plain" name from gsl
            
        ## Install libBigWig
        curl -O https://github.com/dpryan79/libBigWig/archive/refs/tags/0.4.6.tar.gz
        tar -xzf 0.4.6.tar.gz
        mv libBigWig-0.4.6 libBigWig
        cd libBigWig
        make prefix=$INSTALL_DIR install
        cd $INSTALL_DIR
        
        # wiggletools itself (note the modified Makefile)
        curl -O https://github.com/Ensembl/WiggleTools/archive/refs/tags/v1.2.1.tar.gz
        tar -xzf v1.2.1.tar.gz
        mv WiggleTools-1.2.1 WiggleTools
        ./R-3.6.1/bin/Rscript ../scripts/fix_makefile.R
        cd WiggleTools
        make prefix=$INSTALL_DIR
        cd $INSTALL_DIR
    
        #  wigToBigWig -----------------------------------------------------------
        
        if [ $(uname -s) == "Linux" ]; then
            curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
            mv wigToBigWig bin/
        elif [ $(uname -s) == "Darwin" ]; then
            curl -O http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/wigToBigWig
            mv wigToBigWig bin/
        else
            #  Attempt to build from source without assuming the OS
            curl -O http://hgdownload.cse.ucsc.edu/admin/exe/userApps.src.tgz
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
elif [ "$1" == "jhpce" ]; then
    echo "User selected set-up at JHPCE. Installing any missing R packages..."
    module load conda_R/4.1
    Rscript scripts/check_R_packages_JHPCE.R
    
    echo "Setting up test files..."
    Rscript scripts/make_test_manifests.R -d $(pwd)
    
    echo "Preparing main script..."
    sed -i "s|ORIG_DIR=.*|ORIG_DIR=$(pwd)|" run_pipeline_jhpce.sh
    
    echo "Done."
else # neither "docker", "local", nor "jhpce" were chosen
    
    echo 'Error: please specify "docker", "local", or "jhpce" and rerun this script.'
    echo '    eg. bash install_software.sh "local"'
    exit 1
    
fi
