#!/bin/bash

## use qrsh -l mem_free=32G,h_vmem=32G,h_fsize=100G
module load nextflow
export _JAVA_OPTIONS="-Xms25g -Xmx30g"

# The annotation files were downloaded from https://github.com/LieberInstitute/RNAsp/
# specifically from https://github.com/LieberInstitute/RNAsp/tree/438d25b44fef422f67865e4a16f941d50d173399
# and then copied to their location:
#
# git clone git@github.com:LieberInstitute/RNAsp.git
# cd RNAsp
# mv Annotation /dcl01/lieber/ajaffe/lab/RNAsp_static/
# mv Genotyping /dcl01/lieber/ajaffe/lab/RNAsp_static/

## hg38 Simple test (no docker, no sge)
nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	--ercc \
	--fullCov \
	--wg_test \
	-with-report /dcl01/lieber/ajaffe/lab/RNAsp_work/execution_reports/JHPCE_mode_test_run.html \
	-with-dag /dcl01/lieber/ajaffe/lab/RNAsp_work/execution_DAGs/JHPCE_mode_test_run.DAG \
	-profile jhpce,quick \
    -resume \
    -w "/dcl01/lieber/ajaffe/lab/RNAsp_work/runs" \
    --annotation "/dcl01/lieber/ajaffe/lab/RNAsp_static/Annotation" \
    --indexing "/dcl01/lieber/ajaffe/lab/RNAsp_static/Annotation" \
    --genotype "/dcl01/lieber/ajaffe/lab/RNAsp_static/Genotyping" \
    --output "/dcl01/lieber/ajaffe/lab/RNAsp_work/results" \
	$@