#!/bin/bash

## Execution using docker to manage required software, and currently
## configured to use local hardware resources.
export _JAVA_OPTIONS="-Xms8g -Xmx10g"
Software/nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	--ercc \
	--fullCov \
	-with-report execution_reports/docker_run.html \
	-with-dag execution_DAGs/docker_run.html \
	-profile docker \
	$@
