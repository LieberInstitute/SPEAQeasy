#!/bin/bash

## * copy this this script in a working directory that has 
##   a samples.manifest file
## * wrk and results directories are going to be created
## * the output will be in ./results/

### === edit the following lines as needed:
EXPNAME='orgn_ds1' ## experiment label (suffix for result files)
STRAND=reverse # "unstranded" , "forward" or "reverse"
STYPE=paired      # "single" or "paired"
RESUME=1          # set this to 1 if you want to resume interrupted run

SPQZ=/opt/sw/spqz
SPLOG=$PWD/SPEAQeasy_output.log

##REF=hg38 # reference genome - can be hg38, hg19, mm10 or rn6
##--- using custom prepared gencode43 annotation:
REF="hg38main_g43nri"
ANN='/opt/sw/spqz/Annotation/custom'
CUSTOM_ANN="--custom_anno $REF"
##--if you want Gencode25 main chr, uncomment the 3 lines below:
#REF="hg38"
#ANN="/opt/sw/spqz/Annotation"
#CUSTOM_ANN=""

## do not change this unless you need a different profile
PROFILE=localk60 # localk60 uses Gencode 25 main, -k 60 parameter for HISAT2
## Gencode 25 choice is overridden by $CUSTOM_ANN

###=== only edit these if you want to change the default output/input directories
OUTDIR=$(pwd -P)
INDIR="$OUTDIR" # samples.manifest must be here

export NXF_JVM_ARGS="-Xms8g -Xmx10g"

mkdir -p "$OUTDIR/wrk"
mkdir -p "$OUTDIR/results"
optResume=""
if [[ $RESUME == 1 || $RESUME == yes ]]; then
 optResume="-resume"
fi
nextflow -bg -Dnxf.pool.type=sync run $SPQZ/main.nf \
    --sample $STYPE \
    --annotation "$ANN" $CUSTOM_ANN \
    --reference "$REF" \
    --strand $STRAND \
    --trim_mode adaptive \
    --use_salmon \
    --input  "$INDIR" \
         -w  "$OUTDIR/wrk" \
    --output "$OUTDIR/results" \
    --experiment $EXPNAME \
    -profile "$PROFILE" $optResume >& $SPLOG

##removed:  -with-report execution_reports/pipeline_report.html \
