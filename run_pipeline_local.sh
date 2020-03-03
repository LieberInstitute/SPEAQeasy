#!/bin/bash

#  Script to run pipeline using local hardware resources

export _JAVA_OPTIONS="-Xms5g -Xmx7g"

Software/nextflow main.nf \
    --small_test \
    --sample "single" \
    --reference "hg38" \
    --strand "unstranded" \
    --ercc \
    -with-report execution_reports/pipeline_report.html \
    -with-dag execution_DAGs/pipeline_DAG.html \
    -profile local
