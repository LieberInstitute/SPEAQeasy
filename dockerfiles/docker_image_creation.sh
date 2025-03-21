#  This script contains all the code used to build the current docker images,
#  for reproducibility.

#  The user's working directory should be [REPO]/dockerfiles/
REPO=$(pwd)

########################################################
#  R with bioconductor packages pre-installed
########################################################

#  Note that this docker image is deprecated in place of
#  "libddocker: bioc_kallisto:3.14", which works for R 4.1.2 + Bioc 3.14, and
#  includes a few additional R packages. The replacement image is produced
#  here: https://github.com/LieberInstitute/BiocMAP rather than through
#  SPEAQeasy.
mkdir $REPO/r_3.6.1_bioc
cp $REPO/scripts/check_R_packages.R $REPO/r_3.6.1_bioc/
#  produce $REPO/r_3.6.1_bioc/Dockerfile
docker build -t libddocker/r_3.6.1_bioc $REPO/r_3.6.1_bioc
docker push libddocker/r_3.6.1_bioc

########################################################
#  R + bioconductor packages + kallisto
########################################################

#  Note that this docker image is deprecated in place of
#  "libddocker: bioc_kallisto:3.14", which works for R 4.1.2 + Bioc 3.14, and
#  includes a few additional R packages. The replacement image is produced
#  here: https://github.com/LieberInstitute/BiocMAP rather than through
#  SPEAQeasy.
mkdir $REPO/infer_strandness
#  produce $REPO/infer_strandness/Dockerfile
docker build -t libddocker/infer_strandness $REPO/infer_strandness
docker push libddocker/infer_strandness

########################################################
#  regtools
########################################################

#  0.5.1 -----------------------------------------------

#  Originally support python 2.7.9
mkdir -p $REPO/regtools/0.5.1
cp $REPO/scripts/bed_to_juncs.py $REPO/regtools/0.5.1/
# produce $REPO/regtools/0.5.1/Dockerfile
docker build -t libddocker/regtools:0.5.1 $REPO/regtools/0.5.1
docker push libddocker/regtools:0.5.1

#  Update the existing image to only have python 3 and
#  use an updated bed_to_juncs.py
docker build -t libddocker/regtools:0.5.1 $REPO/regtools/0.5.1
docker push libddocker/regtools:0.5.1

#  0.5.33g ---------------------------------------------

mkdir -p $REPO/regtools/0.5.33g
# produce $REPO/regtools/0.5.33g/Dockerfile
docker build -t libddocker/regtools:0.5.33g $REPO/regtools/0.5.33g
docker push libddocker/regtools:0.5.33g

########################################################
#  samtools
########################################################

#  1.9
mkdir -p $REPO/samtools/1.9
#  produce $REPO/samtools/1.9/Dockerfile
docker build -t libddocker/samtools:1.9 $REPO/samtools/1.9
docker push libddocker/samtools:1.9

#  1.10
mkdir -p $REPO/samtools/1.10
#  produce $REPO/samtools/1.10/Dockerfile
docker build -t libddocker/samtools:1.10 $REPO/samtools/1.10
docker push libddocker/samtools:1.9

########################################################
#  samtools + bcftools for variant calling
########################################################

#  1.9
mkdir -p $REPO/variant_calling/1.9
#  produce $REPO/variant_calling/1.9/Dockerfile
docker build -t libddocker/variant_calling:1.9 $REPO/variant_calling/1.9
docker push libddocker/variant_calling:1.9

#  1.10.2
mkdir -p $REPO/variant_calling/1.10.2
#  produce $REPO/variant_calling/1.10.2/Dockerfile
docker build -t libddocker/variant_calling:1.10.2 $REPO/variant_calling/1.10.2
docker push libddocker/variant_calling:1.10.2

########################################################
#  trimmomatic (0.39)
########################################################

mkdir -p $REPO/trimmomatic/0.39
#  produce $REPO/trimmomatic/0.39/Dockerfile
docker build -t libddocker/trimmomatic:0.39 $REPO/trimmomatic/0.39
docker push libddocker/trimmomatic:0.39

########################################################
#  salmon
########################################################

#  0.14.1
mkdir -p $REPO/salmon/0.14.1
#  produce $REPO/salmon/0.14.1/Dockerfile
docker build -t libddocker/salmon:0.14.1 $REPO/salmon/0.14.1
docker push libddocker/salmon:0.14.1

#  Also fixed code in salmon 0.9.1 dockerfile:
cd dockerfiles
docker build -f salmon-0.9.1.dockerfile -t libddocker/salmon-0.9.1:1_v4 .
docker push libddocker/salmon-0.9.1:1_v4
cd ..

#  1.2.1
mkdir -p $REPO/salmon/1.2.1
#  produce $REPO/salmon/1.2.1/Dockerfile
docker build -t libddocker/salmon:1.2.1 $REPO/salmon/1.2.1
docker push libddocker/salmon:1.2.1

########################################################
#  hisat2
########################################################

#  2.1.0
mkdir -p $REPO/hisat2/2.1.0
#  produce $REPO/hisat2/2.1.0/Dockerfile
docker build -t libddocker/hisat2:2.1.0 $REPO/hisat2/2.1.0
docker push libddocker/hisat2:2.1.0

#  2.2.0
mkdir -p $REPO/hisat2/2.2.0
#  produce $REPO/hisat2/2.2.0/Dockerfile
docker build -t libddocker/hisat2:2.2.0 $REPO/hisat2/2.2.0
docker push libddocker/hisat2:2.2.0

#  2.2.1
mkdir -p $REPO/hisat2/2.2.1
#  produce $REPO/hisat2/2.2.1/Dockerfile
docker build -t libddocker/hisat2:2.2.1 $REPO/hisat2/2.2.1
docker push libddocker/hisat2:2.2.1

#  Modify $REPO/hisat2/2.2.1/Dockerfile to include samtools 1.10, in line with
#  the changes making the output from alignment a BAM file rather than SAM
docker build -t libddocker/hisat2:2.2.1 $REPO/hisat2/2.2.1
docker push libddocker/hisat2:2.2.1

########################################################
#  fastQC (0.11.8)
########################################################

mkdir -p $REPO/fastqc/0.11.8
#  produce $REPO/fastqc/0.11.8/Dockerfile
docker build -t libddocker/fastqc:0.11.8 $REPO/fastqc/0.11.8
docker push libddocker/fastqc:0.11.8

########################################################
#  Add "bc" to wiggletools 1.2 image
########################################################

# Add line 7: "RUN apt-get install bc" to dockerfile
cd dockerfiles
docker build -f wiggletools-1.2.dockerfile -t libddocker/wiggletools-1.2:1_v4 .
docker push libddocker/wiggletools-1.2:1_v4
cd ..

########################################################
#  Subread (2.0.0)
########################################################

mkdir -p $REPO/subread/2.0.0
#  produce $REPO/subread/2.0.0/Dockerfile
docker build -t libddocker/subread:2.0.0 $REPO/subread/2.0.0
docker push libddocker/subread:2.0.0

########################################################
#  Kallisto (0.46.1)
########################################################

mkdir -p $REPO/kallisto/0.46.1
#  produce $REPO/kallisto/0.46.1/Dockerfile
docker build -t libddocker/kallisto:0.46.1 $REPO/kallisto/0.46.1
docker push libddocker/kallisto:0.46.1

########################################################
#  RSeQC (3.0.1)
########################################################

mkdir -p $REPO/rseqc/3.0.1
#  produce $REPO/rseqc/3.0.1/Dockerfile
docker build -t libddocker/rseqc:3.0.1 $REPO/rseqc/3.0.1
docker push libddocker/rseqc:3.0.1

########################################################
#  STAR (2.7.8a)
########################################################

mkdir -p $REPO/star/2.7.8a
#  produce $REPO/star/2.7.8a/Dockerfile
docker build -t libddocker/star:2.7.8a $REPO/star/2.7.8a
docker push libddocker/star:2.7.8a
