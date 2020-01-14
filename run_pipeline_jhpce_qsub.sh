#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=50G
#$ -o ./run_pipeline_jhpce_qsub.log
#$ -e ./run_pipeline_jhpce_qsub.log
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow /dcl01/lieber/ajaffe/Nick/RNAsp/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "unstranded" \
    --ercc \
    --small_test \
    -with-report execution_reports/JHPCE_run.html \
    -with-dag execution_DAGs/JHPCE_run.html \
    -profile jhpce \
    --annotation "/users/neagles/rna_sp/Annotation" \
    -w "/dcl01/lieber/ajaffe/lab/RNAsp_work/runs" \
    --output "/dcl01/lieber/ajaffe/lab/RNAsp_work/results" \
	$@