#!/bin/bash

## hg38 Simple test (no docker, no sge)
nextflow main.nf \
	--small_test \
	--sample "single" \
	--reference "hg38" \
	--strand "unstranded" \
	-profile sge,quick \
	-with-report \
	-with-dag flowchart.png \
	$@
