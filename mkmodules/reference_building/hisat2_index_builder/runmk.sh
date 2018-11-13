#!/bin/bash

## This script requests the target to be created by mk
## for this module, the target will be a virtual recipe so we will just name it "HISAT2_INDEX"
## this means that runmk.sh will simply print the order "HISAT2_INDEX" and pass it to mk via xargs
## the command xargs mk $@ ensures that mk can recieve arguments or variables via $@
## for example: bash runmk.sh VARIABLE="MY_VALUE"
echo "HISAT2_INDEX" \
| xargs mk $@
