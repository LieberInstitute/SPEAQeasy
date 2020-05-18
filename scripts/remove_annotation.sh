#!/bin/bash
#$ -cwd

#  Run this script from within the main repository directory! E.g.
#    bash scripts/remove_annotation.sh
#
#  Removes all versions of the specified annotation files, for a
#  given species.
#
#    hisat:       remove assembly fastas and hisat indices
#    transcripts: remove transcript fastas
#    gtf:         remove transcript annotation gtfs

###############################################
#  User configuration
###############################################

species="mm10"

hisat=true
transcripts=false
gtf=false

###############################################

#  Ensure we are in the main repository directory
ls Annotation
if [ $? -ne 0 ]; then
    echo "Exiting- make sure you are in the correct directory when executing this script."
    exit 1
fi

#  Assembly fasta and hisat indices built from it
if [ $hisat == true ]; then
    rm -rf Annotation/reference/$species/assembly/*
fi

#  Transcript annotation gtf
if [ $gtf == true ]; then
    rm -f Annotation/RSeQC/$species/gtf/*
fi

#  Transcript fasta
if [ $transcripts == true ]; then
    rm -rf Annotation/reference/$species/transcripts/*
fi

#  If any of the above are deleted, it implies the annotation objects must be
#  rebuilt (and thus the current ones should be deleted)

if [ $hisat == true ] || [ $gtf == true ] || [ $transcripts == true ]; then
    rm -f Annotation/junction_txdb/chrom_sizes_${species}_*
    rm -f Annotation/junction_txdb/feature_to_Tx_${species}_*.rda
    rm -f Annotation/junction_txdb/junction_annotation_${species}_*.rda
fi
