#!/bin/bash

#bash clean_workdir.sh

## Simple test (no docker, no sge)
#nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -with-report -with-dag flowchart.png

## hg19 Simple test with SGE, (no docker)
#nextflow main.nf --small_test --sample "single" --reference "hg19" --strand "unstranded" --ercc --fullCov -profile sge -with-report -with-dag flowchart.png -resume

## hg38 Simple test with SGE (no docker)
#nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc -profile sge -resume -with-report -with-dag flowchart.png -N iaguilaror@gmail.com

## mm10 Simple test with SGE, (no docker)
#nextflow main.nf --small_test --sample "single" --reference "mm10" --strand "unstranded" --ercc --fullCov -profile sge -with-report -with-dag flowchart.png -resume

## rn6 Simple test with SGE, (no docker)
#nextflow main.nf --small_test --sample "single" --reference "rn6" --strand "unstranded" --ercc --fullCov -profile sge -with-report -with-dag flowchart.png -resume

### Testing reference block of NF pipeline
#nextflow main.nf --sample "single" --reference "hg38" --strand "unstranded" -with-report -with-dag \
#&& nextflow main.nf --sample "single" --reference "hg19" --strand "unstranded" -with-report -with-dag \
#&& nextflow main.nf --sample "single" --reference "mm10" --strand "unstranded" -with-report -with-dag \
#&& nextflow main.nf --sample "single" --reference "rn6" --strand "unstranded" -with-report -with-dag

###Testing hg38 paired unstranded
#nextflow main.nf --small_test --sample "paired" --reference "hg38" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing hg19 paired stranded
#nextflow main.nf --small_test --sample "paired" --reference "hg19" --strand "reverse" --ercc -profile quick -with-report -with-dag flowchart.png
### Testing hg19 paired unstranded
#nextflow main.nf --small_test --sample "paired" --reference "hg19" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing hg19 single stranded
#nextflow main.nf --small_test --sample "single" --reference "hg19" --strand "forward" --ercc -profile quick -with-report -with-dag flowchart.png
### Testing hg19 single unstranded
#nextflow main.nf --small_test --sample "single" --reference "hg19" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png


### Testing hg19 paired stranded
#nextflow main.nf --small_test --sample "single" --reference "hg19" --strand "forward" --ercc -profile quick -with-report -with-dag flowchart.png
### Testing mm10 single unstranded
#nextflow main.nf --small_test --sample "single" --reference "mm10" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png


### Testing mm10 paired unstranded
#nextflow main.nf --small_test --sample "paired" --reference "mm10" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing mm10 paired stranded
#nextflow main.nf --small_test --sample "paired" --reference "mm10" --strand "forward" --ercc -profile quick -with-report -with-dag flowchart.png


### Testing rna6 single unstranded
#nextflow main.nf --small_test --sample "single" --reference "rn6" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing rna6 single stranded
#nextflow main.nf --small_test --sample "single" --reference "rn6" --strand "forward" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing rna6 single stranded
#nextflow main.nf --small_test --sample "single" --reference "rn6" --strand "reverse" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing rna6 paired stranded
#nextflow main.nf --small_test --sample "paired" --reference "rn6" --strand "forward" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing rna6 paired stranded
#nextflow main.nf --small_test --sample "paired" --reference "rn6" --strand "reverse" --ercc -profile quick -with-report -with-dag flowchart.png

### Testing rna6 paired unstranded
nextflow main.nf --small_test --sample "paired" --reference "rn6" --strand "unstranded" --ercc -profile quick -with-report -with-dag flowchart.png
