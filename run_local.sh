#!/bin/bash
## running the pipeline locally on a multi-core machine
## * copy this this script in a working directory that has
##   a samples.manifest file
## Defaults:
## * ./wrk/ and ./results/ are going to be created in the working directory
## * final output will be in ./results/

### === edit the following lines as needed:
STRAND=reverse # "unstranded" , "forward" or "reverse"
STYPE=paired     # "single" or "paired"
RESUME=1          # change this to 1 if you want to resume interrupted run
BG=1    # change this to 1 if you want to run the pipeline in the background
REF=hg38 # reference genome - can be hg38, hg19, mm10 or rn6
### ===

OUTDIR=$(pwd -P)
INDIR="$OUTDIR" # samples.manifest must be here
WRKDIR="$OUTDIR/wrk"

#  After running 'install_local.sh', this should point to the directory
#  where SPEAQeasy was installed and not say "$PWD"
ORIG_DIR=$PWD

#make sure output directories exist
mkdir -p "$WRKDIR"
mkdir -p "$OUTDIR/results"
optResume=""
if [[ $RESUME == 1 || $RESUME == yes ]]; then
 optResume="-resume"
fi
if [[ $BG == 1 || $BG == yes ]]; then
 optBG="-bg"
fi

SPQZ=$ORIG_DIR

SW=$ORIG_DIR/Software
SPLOG=$PWD/SPEAQeasy_output.log

export _JAVA_OPTIONS="-Xms8g -Xmx10g"
export PATH=$SW:$SW/bin:$SW/Trimmomatic-0.39:$SW/FastQC:$PATH

$SPQZ/Software/nextflow -Dnxf.pool.type=sync run $SPQZ/main.nf \
    --sample $STYPE \
    --reference $REF \
    --strand $STRAND \
    --trim_mode adaptive \
    --annotation "$SPQZ/Annotation" \
    --input  "$INDIR" \
         -w  "$WRKDIR" \
    --output "$OUTDIR/results" \
    -with-report execution_reports/pipeline_report.html \
    -profile local $optResume $optBG \
    2>&1 | tee -a $SPLOG


#--- original code below:
#$ORIG_DIR/Software/nextflow $ORIG_DIR/main.nf \
#    --small_test \
#    --sample "single" \
#    --reference "hg38" \
#    --strand "unstranded" \
#    --annotation "$ORIG_DIR/Annotation" \
#    -with-report execution_reports/pipeline_report.html \
#    -profile local \
#    2>&1 | tee -a $SPLOG
    
#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging). Experimental, and may not work on
#  all Linux distributions. Other operating systems are currently not supported.
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
#echo "Generating per-sample logs for debugging..."
#bash $ORIG_DIR/scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
