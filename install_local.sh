#!/bin/bash
#  Usage:  bash install_local.sh
#
# Installation of programs for local running of the pipeline
# Only programs that are not found in the PATH are compiled and installed
#
# Note: versions are not checked !
#
# Assumes core software is installed: java, python3, R, bc, cmake
# Most Linux distros provide packages for these, better use them instead of building them from scratch

errexit() {
    echo -e "Error: $1"
    exit 1
}

trap "errexit" ERR

# must be in the main SPEAQeasy directory
if [[ ! $(fgrep SPEAQeasy main.nf) ]]; then
  errexit "current directory must be the SPEAQeasy code tree !"
fi 

#check core software:
coreProgs=(java python3 bc cmake R Rscript)
for prog in "${coreProgs[@]}" ; do
  if [[ ! -x "$(command -v $prog)" ]]; then
    errexit " core program $prog could not be found.\nPlease make sure it is installed."
  fi
done

SEdir=$(pwd -P)
INSTALL_DIR=$SEdir/Software
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR
export PATH=${INSTALL_DIR}:${INSTALL_DIR}/bin:${INSTALL_DIR}/FastQC:$PATH
#  Install nextflow (latest)
if [[ ! -x "$(command -v nextflow)" ]]; then
  echo "Installing nextflow .. "
  curl -s https://get.nextflow.io | bash
  echo "Please make sure you move the $INSTALL_DIR/nextflow program to a directory in your PATH !"
fi

#  bcftools (1.10.2)  -------------------------------------------------------------
if [[ ! -x "$(command -v bcftools)" ]]; then
 echo "installing bcftools .."
 curl -ksLO https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2
 tar -xjf bcftools.tar.bz2
 cd bcftools-1.15
 ./configure prefix=$INSTALL_DIR
 make
 make install
 cd $INSTALL_DIR
fi

#  fastqc (0.11.8)  -------------------------------------------------------------
if [[ ! -x "$(command -v fastqc)" ]]; then
  echo "installing FastQC .."
  curl -ksLO https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
  unzip fastqc_v0.11.9.zip
  chmod -R 775 FastQC
fi

#  hisat2 (2.2.1)  -------------------------------------------------------------
if [[ ! -x "$(command -v hisat2)" ]]; then
  echo "installing hisat2 .."
  curl -ksLO https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz
  tar -xzf v2.2.1.tar.gz
  cd hisat2-2.2.1
  make -j4
  cd $INSTALL_DIR
  ln -s $INSTALL_DIR/hisat2-2.2.1/hisat2 $INSTALL_DIR/hisat2-2.2.1/hisat2-* .
fi

if [[ ! -x "$(command -v wiggletools)" ]]; then
  echo "installing wiggletools .."
  git clone https://github.com/ebiggers/libdeflate
  cd libdeflate
  make -j 2 libdeflate.a
  cp libdeflate.a ../lib/
  cp libdeflate.h ../include/
  cd $INSTALL_DIR
  curl -ksLO https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2
  tar -xjf htslib-1.15.tar.bz2
  cd htslib-1.15
  ./configure prefix=$INSTALL_DIR
  make
  make install
  cd $INSTALL_DIR
  mv htslib-1.15 htslib # wiggletools expects a "plain" name from htslib
  # install libBigWig
  git clone https://github.com/dpryan79/libBigWig.git
  cd libBigWig
  make prefix=$INSTALL_DIR install
  cd $INSTALL_DIR
  ## Install gsl
  curl -ksLO ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz 
  tar -xvzpf gsl-latest.tar.gz
  /bin/rm -f gsl-latest.tar.gz
  cd gsl*
  ./configure --prefix=$INSTALL_DIR
  make -j4 
  make install
  cd $INSTALL_DIR
  # wiggletools itself (note the modified Makefile)
  curl -o WiggleTools-1.2.11.tar.gz -ksLO https://github.com/Ensembl/WiggleTools/archive/refs/tags/v1.2.11.tar.gz
  tar -xzf WiggleTools-1.2.11.tar.gz
  cd WiggleTools-1.2.11
  #Rscript ../scripts/fix_makefile.R
  perl -i -pe 's/-lBigWig/-l:libBigWig.a/;s/-lhts/-l:libhts.a/;s/-lgsl(\S*)/-l:libgsl$1.a/g;s/-lbz2/-lbz2 -lcrypto -ldeflate/' \
   src/Makefile
  LDFLAGS=-L$INSTALL_DIR/lib CPPFLAGS=-I$INSTALL_DIR/include make prefix=$INSTALL_DIR Wiggletools
  cp bin/wiggletools ..
  cd $INSTALL_DIR
fi

if [[ ! -x "$(command -v samtools)" ]]; then
  curl -ksLO samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
  tar -xjf samtools.tar.bz2
  cd samtools-1.15
  ./configure prefix=$INSTALL_DIR
  make && \
  make install
  cd $INSTALL_DIR
fi

if [[ ! -x "$(command -v kallisto)" ]]; then
  echo "installing kallisto .."
  #  HDF5 (1.10.5), needed for build of kallisto -------------------------------
  curl -O https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
  tar -xzf hdf5-1.10.5.tar.gz
  cd hdf5-1.10.5
  ./configure --disable-parallel --without-szlib --without-pthread --prefix=$INSTALL_DIR
  make
  make install
  cd $INSTALL_DIR
  curl -O https://github.com/pachterlab/kallisto/archive/v0.46.1.tar.gz
  tar -xzf v0.46.1.tar.gz
  cd kallisto-0.46.1
  mkdir build
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
  make
  make prefix=$INSTALL_DIR install
  cd $INSTALL_DIR
fi

#  Install packages that will be used by the pipeline
Rscript ../scripts/check_R_packages.R

BASE_DIR=$(dirname $INSTALL_DIR)
cd $BASE_DIR

#  Signal to load ordinary R packages with 'checkpoint' in each R script
#sed -i "1i #  Added during installation\nlibrary('checkpoint')\ncheckpoint('2021-10-01',\n    project_dir = '$BASE_DIR/scripts/r_packages',\n    checkpoint_location = '$BASE_DIR/Software'\n)\n" scripts/*.R

#  Create the test samples.manifest files
Rscript scripts/make_test_manifests.R -d $(pwd)
cd $INSTALL_DIR

#  regtools (gpertea fork: 0.5.33g)  ----------------------------------
if [[ ! -x "$(command -v regtools)" ]]; then
  curl -O https://github.com/gpertea/regtools/archive/refs/tags/0.5.33g.tar.gz
  tar -xf 0.5.33g.tar.gz
  cd regtools-0.5.33g
  mkdir build
  cd build
  $INSTALL_DIR/bin/cmake ..
  make -j4
  cp regtools $INSTALL_DIR/
  cd $INSTALL_DIR
fi

#  rseqc (3.0.1)  -------------------------------------------------------------

#  Install for the user, because we are not guaranteed system-wide
#  installation privileges
python3 -m pip install --user RSeQC==3.0.1

if [[ 0 -gt 1 ]]; then # not caring about these for now
  #  salmon (1.2.1)  -------------------------------------------------------------
  curl -o salmon_v1.2.1.tar.gz https://github.com/COMBINE-lab/salmon/archive/v1.2.1.tar.gz
  tar -xzf salmon_v1.2.1.tar.gz
  cd salmon-1.2.1
  mkdir build
  cd build
  cmake -DFETCH_BOOST=TRUE ..
  make
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
fi

#  This works on Linux systems but needs testing on MacOS, FreeBSD
if [ $(uname -s) == "Linux" ]; then
    makefile_name="Makefile.Linux"
elif [ $(uname -s) == "Darwin" ]; then
    makefile_name="Makefile.MacOS"
else
    makefile_name="Makefile.FreeBSD"
fi
if [[ ! -x "$(command -v featureCounts)" ]]; then
  curl -ksLO https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-source.tar.gz/download
  tar -zxf download
  cd subread-2.0.3-source/src
  make -f $makefile_name
  cp ../bin/* $INSTALL_DIR/
  cd $INSTALL_DIR
fi

#  trimmomatic (0.39)  -------------------------------------------------------------
if [[ ! -d Trimmomatic-0.39 ]]; then
  curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
  unzip Trimmomatic-0.39.zip
  chmod -R 775 Trimmomatic-0.39
fi

#  wigToBigWig -----------------------------------------------------------
if [[ ! -x "$(command -v wigToBigWig)" ]]; then
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
fi

#  Clean up compressed files
/bin/rm -f $INSTALL_DIR/*.tar.gz
/bin/rm -f $INSTALL_DIR/*.{tgz,bz2,zip}
/bin/rm -f download

#  Fix any strict permissions which would not allow sharing software
#  (and therefore the pipeline as a whole) with those in a group
chmod 775 -R $INSTALL_DIR

#  Point to the original repository so that the "main" scripts can be
#  trivially copied to share the pipeline
sed -i "s|ORIG_DIR=\$PWD|ORIG_DIR=$SEdir|" ../run_local.sh
sed -i "s|ORIG_DIR=\$PWD|ORIG_DIR=$SEdir|" ../run_pipeline_sge.sh
sed -i "s|ORIG_DIR=\$PWD|ORIG_DIR=$SEdir|" ../run_pipeline_local.sh
sed -i "s|ORIG_DIR=\$PWD|ORIG_DIR=$SEdir|" ../run_pipeline_slurm.sh
