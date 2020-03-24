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
    -profile local \
    > SPEAQeasy_output.log
    
#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging). Experimental, and may not work on
#  all Linux distributions. Other operating systems are currently not supported.
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
