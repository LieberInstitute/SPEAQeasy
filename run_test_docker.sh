#!/bin/bash

## hg38 Simple test (no docker, no sge)
nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	--ercc \
	--fullCov \
        -profile docker,quick \
	-with-report \
	-with-dag flowchart.png \
	-resume \
	$1
