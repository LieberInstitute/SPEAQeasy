#!/bin/bash

#bash clean_workdir.sh

## Simple test (no docker, no sge)
#nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov

## Simple test with SGE (no docker)
nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -profile sge
