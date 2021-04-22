# SPEAQeasy- an RNA-seq Pipeline from the Jaffe Lab at Lieber Institute #

## Summary ##

SPEAQeasy is a **S**calable RNA-seq **P**ipeline for **E**xpression **A**nalysis and **Q**uantification based on the [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline). Built on [nextflow](https://www.nextflow.io/), and capable of using Docker containers and utilizing common resource managers (e.g. SLURM), this port of the RNAseq-pipeline can be used in different computer environments. It is described in the following bioRxiv pre-print: https://www.biorxiv.org/content/10.1101/2020.12.11.386789v1.

The main function of this pipeline is to produce comparable files to those used in [recount2](https://jhubiostatistics.shinyapps.io/recount/), a tool that provides gene, exon, exon-exon junction and base-pair level data.

This pipeline allows researchers to contribute data to the recount2 project even from outside the [JHPCE](https://jhpce.jhu.edu/).


## Workflow overview ##

![General Workflow](https://github.com/LieberInstitute/SPEAQeasy/blob/master/notes/workflow.png)

SPEAQeasy takes raw RNA-seq reads and produces analysis-ready R objects, providing a "bridge to the Bioconductor universe", where researchers can utilize the powerful existing set of tools to quickly perform desired analyses.

Beginning with a set of FASTQ files (optionally gzipped), SPEAQeasy ultimately produces [`RangedSummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) objects to store gene, exon, and exon-exon junction counts for an experiment. Optionally, expressed regions data is generated, enabling easy computation of differentially expressed regions (DERs).

[Our vignette](http://research.libd.org/SPEAQeasy-example) demonstrates how genotype calls by SPEAQeasy can be coupled with user-provided genotype and phenotype data to easily resolve identity issues that arise during sequencing. We then walk through an example differential expression analysis and explore data visualization options.

## Pipeline features ##

- Automatically merge samples split across multiple FASTQ files, using the `samples.manifest` [input](http://research.libd.org/SPEAQeasy/manifest.html)
- Trivially [select any GENCODE annotation release](http://research.libd.org/SPEAQeasy/annotation.html#choosing-a-release) for "hg38", "hg19", or "mm10" references (Ensembl for "rn6" reference) and [adjust other annotation settings](http://research.libd.org/SPEAQeasy/annotation.html#choosing-a-build) with simple configuration
- Generates a single VCF file for experiments on human reference, which can be used to [resolve sample identity issues](http://research.libd.org/SPEAQeasy-example) and salvage problematic samples
- Supports [docker to manage software dependencies and is preconfigured for execution locally or on SLURM or SGE clusters](http://research.libd.org/SPEAQeasy/setup-details.html#installation)
- Multiple users can [share a single SPEAQeasy installation](http://research.libd.org/SPEAQeasy/setup-details.html#sharing) with minimal work
- Detailed, [user-friendly logging](http://research.libd.org/SPEAQeasy/help.html#checking-per-sample-logs) for transparency and identifying potential issues


## Getting started ##

The [SPEAQeasy documentation website](http://research.libd.org/SPEAQeasy/index.html) describes the pipeline in full detail. For briefly getting started, check out the [quick start guide](http://research.libd.org/SPEAQeasy/quick-start.html).

Because SPEAQeasy is based on the [nextflow](https://www.nextflow.io/) workflow manager, it supports execution on computing clusters managed by [SLURM](https://slurm.schedmd.com/overview.html) or [SGE](https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html) without any configuration (local execution is also possible). Those with access to [docker](https://www.docker.com/) can very simply use docker containers to manage SPEAQeasy software dependencies, though we provide a script for installing dependencies for users without docker or even root privileges.


## Authors ##

Original Pipeline

 [Emily Burke](mailto:emily.burke@libd.org>),
 [Leonardo Collado-Tores](mailto:lcolladotor@gmail.com),
 [Andrew Jaffe](mailto:andrew.jaffe@libd.org),
 [BaDoi Phan](mailto:badoi.phan@pitt.edu) 
 
Nextflow Port

 [Nick Eagles](mailto:nick.eagles@libd.org),
 [Brianna Barry](https://github.com/BriannaBarry),
 [Jacob Leonard](mailto:leonard.jacob09@gmail.com),
 [Israel Aguilar](mailto:iaguilaror@gmail.com),
 [Violeta Larios](mailto:siedracko@gmail.com),
 [Everardo Gutierrez](mailto:ever.gmillan@gmail.com)

## Contact ##

* [Leonardo Collado Torres](http://lcolladotor.github.io/)
* [Winter Genomics Team](http://www.wintergenomics.com)
