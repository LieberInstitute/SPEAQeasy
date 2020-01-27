#  This script contains all the code used to build the current docker images,
#  for reproducibility.

#  The user's working directory should be [REPO]/dockerfiles/
REPO=..

########################################################
#  R with bioconductor packages pre-installed
########################################################

mkdir $REPO/r_3.6.1_bioc
cp $REPO/scripts/check_R_packages.R $REPO/r_3.6.1_bioc/
#  produce $REPO/r_3.6.1_bioc/Dockerfile
docker build -t libddocker/r_3.6.1_bioc $REPO/r_3.6.1_bioc
docker push libddocker/r_3.6.1_bioc

########################################################
#  R + bioconductor packages + kallisto
########################################################

mkdir $REPO/infer_strandness
#  produce $REPO/infer_strandness/Dockerfile
docker build -t libddocker/infer_strandness $REPO/infer_strandness
docker push libddocker/infer_strandness

########################################################
#  regtools (0.5.1) and python 2.7.9
########################################################

mkdir -p $REPO/regtools/0.5.1
cp $REPO/scripts/bed_to_juncs.py $REPO/regtools/0.5.1/
# produce $REPO/regtools/0.5.1/Dockerfile
docker build -t libddocker/regtools:0.5.1 $REPO/regtools/0.5.1
docker push libddocker/regtools:0.5.1

########################################################
#  samtools (1.9)
########################################################

mkdir -p $REPO/samtools/1.9
#  produce $REPO/samtools/1.9/Dockerfile
docker build -t libddocker/samtools:1.9 $REPO/samtools/1.9
docker push libddocker/samtools:1.9

########################################################
#  samtools + bcftools for variant calling
########################################################

mkdir -p $REPO/variant_calling/1.9
#  produce $REPO/variant_calling/1.9/Dockerfile
docker build -t libddocker/variant_calling:1.9 $REPO/variant_calling/1.9
docker push libddocker/variant_calling:1.9