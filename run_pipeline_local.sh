#!/bin/bash

## Test using local hardware resources and software pointed to
## in conf/command_paths_long.config
export _JAVA_OPTIONS="-Xms5g -Xmx7g"
./Software/nextflow main.nf \
    --small_test \
    --sample "single" \
    --reference "hg38" \
    --strand "unstranded" \
    --ercc \
    --fullCov \
    --force_trim \
    -with-report execution_reports/System_mode_test_run.html \
    -with-dag execution_DAGs/System_mode_test_run.png \
    -profile standard \
    $@
