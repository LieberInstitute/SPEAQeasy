#!/bin/bash

#bash clean_workdir.sh

## Simple test (no docker, no sge)
#nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -with-report -with-dag flowchart.png

## Simple test with SGE (no docker)
#nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -profile sge -resume -with-report -with-dag flowchart.png

## Simple test (no docker, no sge) Single_stranded
#nextflow main.nf --sample "single" --reference "hg38" --strand "reverse" --ercc --fullCov -with-report -with-dag flowchart.png

## Simple test (no docker, no sge) Paired_stranded forward
#nextflow main.nf --sample "paired" --reference "hg38" --strand "forward" --ercc --fullCov -with-report -with-dag flowchart.png

## Simple test (no docker, no sge) Paired_stranded reverse
#nextflow main.nf --sample "paired" --reference "hg38" --strand "reverse" --ercc --fullCov -with-report -with-dag flowchart.png

## Simple test (no docker, no sge) Paired_stranded reverse
nextflow main.nf --sample "paired" --reference "hg38" --strand "unstranded" --ercc --fullCov -with-report -with-dag flowchart.png -resume

