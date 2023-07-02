#!/bin/bash
#$ -l bluejay,mem_free=28G,h_vmem=32G,h_fsize=800G
#$ -o ./SPEAQeasy_output.log
#$ -e ./SPEAQeasy_output.log
#$ -cwd

## * copy this this script in a working directory that has 
##   a samples.manifest file
## * wrk and results directories are going to be created
## * the output will be in ./results/

### === edit the following lines as needed:
STRAND=reverse # "unstranded" , "forward" or "reverse"
STYPE=paired    # "single" or "paired"
RESUME=0        # change this to 1 if you want to resume interrupted run
REF=hg38 # reference genome - can be hg38, hg19, mm10 or rn6

### === only edit these if you want to change the default output/input directories
OUTDIR=$(pwd -P)
INDIR="$OUTDIR" # samples.manifest must be here

#module use /jhpce/shared/jhpce/modulefiles/libd
module load nextflow/20.01.0
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

SPQZ=/dcl02/lieber/ajaffe/gpertea/spqz
#SPQZ=/dcs04/lieber/lcolladotor/dbDev_LIBD001/RNAseq/SPEAQeasy
## /dcl01/lieber/ajaffe/Nick/SPEAQeasy
SPLOG=$PWD/SPEAQeasy_output.log

mkdir -p "$OUTDIR/wrk"
mkdir -p "$OUTDIR/results"
optResume=""
if [[ $RESUME == 1 || $RESUME == yes ]]; then
 optResume="-resume"
fi
## use the Annotation directly from Nick's installation?
## ANN="/dcl01/lieber/ajaffe/Nick/SPEAQeasy/Annotation"
ANN="$SPQZ/Annotation"
## for local: nextflow -bg -Dnxf.pool.type=sync run $SPQZ/main.nf
nextflow $SPQZ/main.nf \
    --sample $STYPE \
    --reference $REF \
    --strand $STRAND \
    --trim_mode adaptive \
    --annotation "$ANN" \
    --input  "$INDIR" \
         -w  "$OUTDIR/wrk" \
    --output "$OUTDIR/results" \
    -with-report execution_reports/pipeline_report.html \
    -profile jhpce_G25pri $optResume
## for local: -profile local $optResume
