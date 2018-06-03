#!/bin/bash

## hg38 Simple test (no docker, no sge)
nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	--ercc \
	--fullCov \
	--wg_test \
	-with-report execution_reports/System_mode_test_run.html \
	-with-dag execution_DAGs/System_mode_test_run.png \
	-resume \
	-profile docker,quick \
	$@
