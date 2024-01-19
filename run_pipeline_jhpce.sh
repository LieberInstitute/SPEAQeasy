#!/bin/bash

#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH --job-name=SPEAQeasy
#SBATCH -o ./SPEAQeasy_output.log
#SBATCH -e ./SPEAQeasy_output.log

#  After running 'install_software.sh', this should point to the directory
#  where SPEAQeasy was installed, and not say "$PWD"
ORIG_DIR=$PWD

module load nextflow/23.10.0
export NXF_JVM_ARGS="-Xms8g -Xmx10g"

nextflow run $ORIG_DIR/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "forward" \
    --small_test \
    --annotation "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation" \
    -with-report execution_reports/JHPCE_run.html \
    -profile jhpce

#   Log successful runs on non-test data in a central location. Please adjust
#   the log path here if it is changed at the top!
bash $ORIG_DIR/scripts/track_runs.sh $PWD/SPEAQeasy_output.log

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $ORIG_DIR/scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
