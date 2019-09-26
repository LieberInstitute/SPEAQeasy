#!/bin/bash
#$ -l mem_free=60G,h_vmem=60G,h_fsize=50G
#$ -o ./run_test_jhpce_qsub_out.txt
#$ -e ./run_test_jhpce_qsub_out.txt
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow main.nf \
	--sample "paired" \
	--reference "hg19" \
	--strand "unstranded" \
	--ercc \
  -resume \
 	-with-report execution_reports/JHPCE_mode_test_run.html \
	-with-dag execution_DAGs/JHPCE_mode_test_run.DAG \
	-profile jhpce \
  -w "/dcl01/lieber/ajaffe/lab/RNAsp_work/runs" \
  --input "/users/neagles/rna_sp/input/" \
  --annotation "/dcl01/lieber/ajaffe/lab/RNAsp_static/Annotation"  \
  --genotype "/dcl01/lieber/ajaffe/lab/RNAsp_static/Genotyping" \
  --output "/dcl01/lieber/ajaffe/lab/RNAsp_work/results" \
	$@