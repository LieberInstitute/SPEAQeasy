#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=800G
#$ -o ./SPEAQeasy_output.log
#$ -e ./SPEAQeasy_output.log
#$ -cwd

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "forward" \
    --small_test \
    --annotation "/dcl01/lieber/ajaffe/Nick/SPEAQeasy/Annotation" \
    -with-report execution_reports/JHPCE_run.html \
    -profile jhpce

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
