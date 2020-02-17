#!/bin/bash

## Script to run the pipeline in a Sun Grid Engines (SGE) environment

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

Software/nextflow main.nf \
    --small_test \
    --sample "single" \
    --reference "hg38" \
    --strand "unstranded" \
    --ercc \
    -with-report execution_reports/pipeline_report.html \
    -with-dag execution_DAGs/pipeline_DAG.html \
    -profile sge
