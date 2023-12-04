#!/bin/bash
#$ -l mem_free=18G,h_vmem=28G,h_fsize=800G
#$ -o ./SPEAQeasy_output.log
#$ -e ./SPEAQeasy_output.log
#$ -cwd

## * copy this this script in a working directory that has 
##   a samples.manifest file
## * wrk and results directories are going to be created
## * the output will be in ./results/

### === edit the following lines as needed:
STRAND=reverse # "unstranded" , "forward" or "reverse"
STYPE=paired      # "single" or "paired"
RESUME=1          # change this to 1 if you want to resume interrupted run

SPQZ=/opt/sw/spqz
SPLOG=$PWD/SPEAQeasy_output.log

##REF=hg38 # reference genome - can be hg38, hg19, mm10 or rn6
##--- using custom gencode43 nri annotation:
REF="hg38main_g43nri"
ANN=/opt/sw/spqz/Annotation/custom

### === only edit these if you want to change the default output/input directories
OUTDIR=$(pwd -P)
INDIR="$OUTDIR" # samples.manifest must be here

#export _JAVA_OPTIONS="-Xms8g -Xmx10g"
export NXF_JVM_ARGS="-Xms8g -Xmx10g"

mkdir -p "$OUTDIR/wrk"
mkdir -p "$OUTDIR/results"
optResume=""
if [[ $RESUME == 1 || $RESUME == yes ]]; then
 optResume="-resume"
fi
nextflow -bg -Dnxf.pool.type=sync run $SPQZ/main.nf \
    --sample $STYPE \
    --annotation "$ANN" \
    --custom_anno "$REF" \
    --reference "$REF" \
    --strand $STRAND \
    --trim_mode adaptive \
    --use_salmon \
    --input  "$INDIR" \
         -w  "$OUTDIR/wrk" \
    --output "$OUTDIR/results" \
    --experiment 'orgn_ds2' \
    -profile localk60 $optResume >& $SPLOG

##removed:  -with-report execution_reports/pipeline_report.html \
