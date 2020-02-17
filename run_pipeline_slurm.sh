#!/bin/bash

## Script to run the pipeline on a SLURM cluster

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/nextflow main.nf \
    --sample "single" \
    --reference "hg19" \
    --strand "unstranded" \
    --small_test \
    -with-report execution_reports/pipeline_report.html \
    -with-dag execution_DAGs/pipeline_DAG.html \
    -profile slurm
