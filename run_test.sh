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
	-with-report test_report.html \
	-with-dag test_flowchart.png \
	-resume \
	-profile quick \
	$@
