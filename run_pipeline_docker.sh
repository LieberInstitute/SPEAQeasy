#!/bin/bash

## Execution using docker to manage required software, and currently
## configured to use local hardware resources.
export _JAVA_OPTIONS="-Xms8g -Xmx10g"
nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	--ercc \
	--fullCov \
	-with-report execution_reports/Docker_mode_test_run.html \
	-with-dag execution_DAGs/Docker_mode_test_run.png \
	-resume \
	-profile docker,quick \
	$@
