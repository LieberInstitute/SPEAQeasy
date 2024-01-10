#!/bin/bash
#$ -l virtual_free=40G,h_fsize=800G
#$ -o SPEAQeasy_output.log
#$ -e SPEAQeasy_output.log
#$ -cwd

## Script to run the pipeline in a Sun Grid Engines (SGE) environment

#  After running 'install_software.sh', this should point to the directory
#  where SPEAQeasy was installed, and not say "$PWD"
ORIG_DIR=$PWD

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

$ORIG_DIR/Software/nextflow run $ORIG_DIR/main.nf \
    --small_test \
    --sample "single" \
    --reference "hg38" \
    --strand "unstranded" \
    --annotation "$ORIG_DIR/Annotation" \
    -with-report execution_reports/pipeline_report.html \
    -profile sge

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging). Experimental, and may not work on
#  all Linux distributions. Other operating systems are currently not supported.
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $ORIG_DIR/scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
