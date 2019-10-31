#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=50G
#$ -o ./run_test_jhpce_qsub_out.txt
#$ -e ./run_test_jhpce_qsub_out.txt
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "unstranded" \
    --ercc \
    --small_test \
    -with-report execution_reports/JHPCE_mode_test_run.html \
    -with-dag execution_DAGs/JHPCE_mode_test_run.DAG \
    -profile jhpce \
    --annotation "/users/neagles/rna_sp/Annotation" \
    -w "/dcl01/lieber/ajaffe/lab/RNAsp_work/runs" \
    --output "/dcl01/lieber/ajaffe/lab/RNAsp_work/results" \
	$@