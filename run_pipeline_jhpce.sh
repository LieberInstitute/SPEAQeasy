#!/bin/bash

## use 'qrsh -l mem_free=60G,h_vmem=60G,h_fsize=100G' for real runs
module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

# The annotation files were downloaded from https://github.com/LieberInstitute/RNAsp/
# specifically from https://github.com/LieberInstitute/RNAsp/tree/438d25b44fef422f67865e4a16f941d50d173399
# and then copied to their location:
#
# git clone git@github.com:LieberInstitute/RNAsp.git
# cd RNAsp
# mv Annotation /dcl01/lieber/ajaffe/lab/RNAsp_static/
# mv Genotyping /dcl01/lieber/ajaffe/lab/RNAsp_static/

nextflow /dcl01/lieber/ajaffe/Nick/RNAsp/main.nf \
    --sample "single" \
    --reference "hg19" \
    --strand "unstranded" \
    --ercc \
    --small_test \
    -with-report /dcl01/lieber/ajaffe/lab/RNAsp_work/execution_reports/JHPCE_run.html \
    -with-dag /dcl01/lieber/ajaffe/lab/RNAsp_work/execution_DAGs/JHPCE_run.html \
    -profile jhpce  \
    -w "/dcl01/lieber/ajaffe/lab/RNAsp_work/runs" \
    --annotation "/users/neagles/rna_sp/Annotation" \
    --output "/dcl01/lieber/ajaffe/lab/RNAsp_work/results" \
    $@

