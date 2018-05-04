#!/bin/bash

bash clean_workdir.sh

nextflow main.nf --sample "single" --reference "hg38" --strand "unstranded" --small_test --merge
