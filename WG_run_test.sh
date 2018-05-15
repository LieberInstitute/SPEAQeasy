#!/bin/bash

#bash clean_workdir.sh

## Simple test (no docker, no sge)
#nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -with-report -with-dag flowchart.png

## hg19 Simple test with SGE, (no docker)
#nextflow main.nf --small_test --sample "single" --reference "hg19" --strand "unstranded" --ercc --fullCov -profile sge -with-report -with-dag flowchart.png -resume

## hg38 Simple test with SGE (no docker)
nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc -profile sge -resume -with-report -with-dag flowchart.png -N iaguilaror@gmail.com

## mm10 Simple test with SGE, (no docker)
#nextflow main.nf --small_test --sample "single" --reference "mm10" --strand "unstranded" --ercc --fullCov -profile sge -with-report -with-dag flowchart.png -resume

## rn6 Simple test with SGE, (no docker)
#nextflow main.nf --small_test --sample "single" --reference "rn6" --strand "unstranded" --ercc --fullCov -profile sge -with-report -with-dag flowchart.png -resume

### Testing reference block of NF pipeline
#nextflow main.nf --sample "single" --reference "hg38" --strand "unstranded" -with-report -with-dag \
#&& nextflow main.nf --sample "single" --reference "hg19" --strand "unstranded" -with-report -with-dag \
#&& nextflow main.nf --sample "single" --reference "mm10" --strand "unstranded" -with-report -with-dag \
#&& nextflow main.nf --sample "single" --reference "rn6" --strand "unstranded" -with-report -with-dag

