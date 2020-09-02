#!/bin/bash
#SBATCH --output=SPEAQeasy_output.log
#SBATCH --mem=40G

## Script to run the pipeline on a SLURM cluster

#  After running 'install_software.sh', this should point to the directory
#  where SPEAQeasy was installed, and not say "$PWD"
ORIG_DIR=$PWD

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

$ORIG_DIR/Software/nextflow main.nf \
    --sample "single" \
    --reference "hg19" \
    --strand "unstranded" \
    --small_test \
    --annotation "$ORIG_DIR/Annotation" \
    -with-report execution_reports/pipeline_report.html \
    -profile slurm
    
#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging). Experimental, and may not work on
#  all Linux distributions. Other operating systems are currently not supported.
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $ORIG_DIR/scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
