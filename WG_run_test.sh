#!/bin/bash

#bash clean_workdir.sh

#################################################
## System mode test (no docker, no sge)
############################################
#### hg38 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg38" --sample "single" --strand "unstranded" \
#-with-report execution_reports/System_mode_noDocker_noSGE_hg38_nomerge_single_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg38_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg38" --sample "single" --strand "forward" \
#-with-report execution_reports/System_mode_noDocker_noSGE_hg38_nomerge_single_forward.html -with-dag execution_DAGs/execution_DAGsSystem_mode_noDocker_noSGE_hg38_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg38" --sample "single" --strand "reverse" \
#-with-report execution_reports/System_mode_noDocker_noSGE_hg38_nomerge_single_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg38_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg38" --sample "paired" --strand "unstranded" \
-with-report execution_reports/System_mode_noDocker_noSGE_hg38_nomerge_paired_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg38_nomerge_paired_unstranded.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg38" --sample "paired" --strand "forward" \
-with-report execution_reports/System_mode_noDocker_noSGE_hg38_nomerge_paired_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg38_nomerge_paired_forward.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg38" --sample "paired" --strand "reverse" \
-with-report execution_reports/System_mode_noDocker_noSGE_hg38_nomerge_paired_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg38_nomerge_paired_reverse.png

#### mm10 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "mm10" --sample "single" --strand "unstranded" \
#-with-report execution_reports/System_mode_noDocker_noSGE_mm10_nomerge_single_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_mm10_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "mm10" --sample "single" --strand "forward" \
#-with-report execution_reports/System_mode_noDocker_noSGE_mm10_nomerge_single_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_mm10_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "mm10" --sample "single" --strand "reverse" \
#-with-report execution_reports/System_mode_noDocker_noSGE_mm10_nomerge_single_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_mm10_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "mm10" --sample "paired" --strand "unstranded" \
-with-report execution_reports/System_mode_noDocker_noSGE_mm10_nomerge_paired_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_mm10_nomerge_paired_unstranded.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "mm10" --sample "paired" --strand "forward" \
-with-report execution_reports/System_mode_noDocker_noSGE_mm10_nomerge_paired_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_mm10_nomerge_paired_forward.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "mm10" --sample "paired" --strand "reverse" \
-with-report execution_reports/System_mode_noDocker_noSGE_mm10_nomerge_paired_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_mm10_nomerge_paired_reverse.png

#### rn6 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "rn6" --sample "single" --strand "unstranded" \
#-with-report execution_reports/System_mode_noDocker_noSGE_rn6_nomerge_single_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_rn6_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "rn6" --sample "single" --strand "forward" \
#-with-report execution_reports/System_mode_noDocker_noSGE_rn6_nomerge_single_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_rn6_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "rn6" --sample "single" --strand "reverse" \
#-with-report execution_reports/System_mode_noDocker_noSGE_rn6_nomerge_single_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_rn6_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "rn6" --sample "paired" --strand "unstranded" \
-with-report execution_reports/System_mode_noDocker_noSGE_rn6_nomerge_paired_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_rn6_nomerge_paired_unstranded.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "rn6" --sample "paired" --strand "forward" \
-with-report execution_reports/System_mode_noDocker_noSGE_rn6_nomerge_paired_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_rn6_nomerge_paired_forward.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "rn6" --sample "paired" --strand "reverse" \
-with-report execution_reports/System_mode_noDocker_noSGE_rn6_nomerge_paired_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_rn6_nomerge_paired_reverse.png

#### hg19 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg19" --sample "single" --strand "unstranded" \
#-with-report execution_reports/System_mode_noDocker_noSGE_hg19_nomerge_single_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg19_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg19" --sample "single" --strand "forward" \
#-with-report execution_reports/System_mode_noDocker_noSGE_hg19_nomerge_single_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg19_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg19" --sample "single" --strand "reverse" \
#-with-report execution_reports/System_mode_noDocker_noSGE_hg19_nomerge_single_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg19_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg19" --sample "paired" --strand "unstranded" \
-with-report execution_reports/System_mode_noDocker_noSGE_hg19_nomerge_paired_unstranded.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg19_nomerge_paired_unstranded.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg19" --sample "paired" --strand "forward" \
-with-report execution_reports/System_mode_noDocker_noSGE_hg19_nomerge_paired_forward.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg19_nomerge_paired_forward.png

nextflow main.nf  -resume -profile quick --small_test --ercc --fullCov --reference "hg19" --sample "paired" --strand "reverse" \
-with-report execution_reports/System_mode_noDocker_noSGE_hg19_nomerge_paired_reverse.html -with-dag execution_DAGs/System_mode_noDocker_noSGE_hg19_nomerge_paired_reverse.png

#################################################
## SGE mode test (no docker)
############################################
#### hg38 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg38" --sample "single" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_hg38_nomerge_single_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_hg38_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg38" --sample "single" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_hg38_nomerge_single_forward.html -with-dag execution_DAGs/execution_DAGsSGE_mode_noDocker_hg38_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg38" --sample "single" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_hg38_nomerge_single_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_hg38_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg38" --sample "paired" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_hg38_nomerge_paired_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_hg38_nomerge_paired_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg38" --sample "paired" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_hg38_nomerge_paired_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_hg38_nomerge_paired_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg38" --sample "paired" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_hg38_nomerge_paired_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_hg38_nomerge_paired_reverse.png

#### mm10 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "mm10" --sample "single" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_mm10_nomerge_single_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_mm10_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "mm10" --sample "single" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_mm10_nomerge_single_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_mm10_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "mm10" --sample "single" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_mm10_nomerge_single_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_mm10_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "mm10" --sample "paired" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_mm10_nomerge_paired_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_mm10_nomerge_paired_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "mm10" --sample "paired" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_mm10_nomerge_paired_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_mm10_nomerge_paired_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "mm10" --sample "paired" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_mm10_nomerge_paired_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_mm10_nomerge_paired_reverse.png

#### rn6 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "rn6" --sample "single" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_rn6_nomerge_single_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_rn6_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "rn6" --sample "single" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_rn6_nomerge_single_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_rn6_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "rn6" --sample "single" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_rn6_nomerge_single_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_rn6_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "rn6" --sample "paired" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_rn6_nomerge_paired_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_rn6_nomerge_paired_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "rn6" --sample "paired" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_rn6_nomerge_paired_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_rn6_nomerge_paired_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "rn6" --sample "paired" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_rn6_nomerge_paired_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_rn6_nomerge_paired_reverse.png

#### hg19 block ####

## singles, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg19" --sample "single" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_hg19_nomerge_single_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_hg19_nomerge_single_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg19" --sample "single" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_hg19_nomerge_single_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_hg19_nomerge_single_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg19" --sample "single" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_hg19_nomerge_single_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_hg19_nomerge_single_reverse.png

## paired, no docker, no SGE, no Merge
#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg19" --sample "paired" --strand "unstranded" \
#-with-report execution_reports/SGE_mode_noDocker_hg19_nomerge_paired_unstranded.html -with-dag execution_DAGs/SGE_mode_noDocker_hg19_nomerge_paired_unstranded.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg19" --sample "paired" --strand "forward" \
#-with-report execution_reports/SGE_mode_noDocker_hg19_nomerge_paired_forward.html -with-dag execution_DAGs/SGE_mode_noDocker_hg19_nomerge_paired_forward.png

#nextflow main.nf  -resume -profile quick,sge --small_test --ercc --fullCov --reference "hg19" --sample "paired" --strand "reverse" \
#-with-report execution_reports/SGE_mode_noDocker_hg19_nomerge_paired_reverse.html -with-dag execution_DAGs/SGE_mode_noDocker_hg19_nomerge_paired_reverse.png
