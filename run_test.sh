#!/bin/bash

bash clean_workdir.sh

nextflow main.nf --small_test --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov
