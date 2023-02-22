# SPEAQeasy- a **S**calable **P**ipeline for **E**xpression **A**nalysis and **Q**uantification that is **easy** to install and share #

## Summary ##

SPEAQeasy is a **S**calable RNA-seq **P**ipeline for **E**xpression **A**nalysis and **Q**uantification based on the [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline). Built on [nextflow](https://www.nextflow.io/), and capable of using Docker containers and utilizing common resource managers (e.g. SLURM), this port of the RNAseq-pipeline can be used in different computer environments. It is described in the manuscript [here](https://doi.org/10.1186/s12859-021-04142-3).

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
 
## Cite `SPEAQeasy`

We hope that [`SPEAQeasy`](http://research.libd.org/SPEAQeasy/) will be useful for your research. Please use the following bibtex information to cite the software and overall approach. Thank you!

```
@article {Eagles2021,
	author = {Eagles, Nicholas J. and Burke, Emily E. and Leonard, Jabob and Barry, Brianna K. and Stolz, Joshua M. and Huuki, Louise and Phan, BaDoi N. and Larrios Serrato, Violeta and Guti{\'e}rrez-Mill{\'a}n, Everardo and Aguilar-Ordo{\~n}ez, Israel and Jaffe, Andrew E. and Collado-Torres, Leonardo},
	title = {SPEAQeasy: a scalable pipeline for expression analysis and quantification for R/bioconductor-powered RNA-seq analyses},
	year = {2021},
	doi = {10.1186/s12859-021-04142-3},
	publisher = {Springer Science and Business Media LLC},
	URL = {https://doi.org/10.1186/s12859-021-04142-3},
	journal = {BMC Bioinformatics}
}
```

## Contact ##

* [Leonardo Collado Torres](http://lcolladotor.github.io/)
* [Winter Genomics Team](http://www.wintergenomics.com)
