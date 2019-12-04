#!/bin/bash

## Test script to run the pipeline in a Sun Grid Engines (SGE) environment
export _JAVA_OPTIONS="-Xms8g -Xmx10g"
nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	--ercc \
	--fullCov \
	-with-report execution_reports/SGE_mode_test_run.html \
	-with-dag execution_DAGs/SGE_mode_test_run.png \
	-resume \
	-profile sge,quick \
	$@
