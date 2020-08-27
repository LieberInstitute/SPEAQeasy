# SPEAQeasy- an RNA-seq Pipeline from the Jaffe Lab at Lieber Institute #

## Summary ##

SPEAQeasy is a **S**calable RNA-seq **P**ipeline for **E**xpression **A**nd **Q**uantification based on the [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline). Built on [nextflow](https://www.nextflow.io/), and capable of using Docker containers and utilizing common resource managers (e.g. SLURM), this port of the RNAseq-pipeline can be used in different computer environments.

The main function of this pipeline is to produce comparable files to those used in [recount2](https://jhubiostatistics.shinyapps.io/recount/), a tool that provides gene, exon, exon-exon junction and base-pair level data.

This pipeline allows researchers to contribute data to the recount2 project even from outside the [JHPCE](https://jhpce.jhu.edu/).


## Workflow overview ##

![General Workflow](https://github.com/LieberInstitute/SPEAQeasy/blob/master/notes/workflow.png)

SPEAQeasy takes raw RNA-seq reads and produces analysis-ready R objects, providing a "bridge to the Bioconductor universe", where researchers can utilize the powerful existing set of tools to quickly perform desired analyses.

Beginning with a set of FASTQ files (optionally gzipped), SPEAQeasy ultimately produces [`RangedSummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) objects to store gene, exon, and exon-exon junction counts for an experiment. Optionally, expressed regions data is generated, enabling easy computation of differentially expressed regions (DERs).

[Our vignette](http://research.libd.org/SPEAQeasy-example) demonstrates how genotype calls by SPEAQeasy can be coupled with user-provided genotype and phenotype data to easily resolve identity issues that arise during sequencing. We then walk through an example differential expression analysis and explore data visualization options.


## Getting started ##

The [SPEAQeasy documentation website](http://research.libd.org/SPEAQeasy/index.html) describes the pipeline in full detail. For briefly getting started, check out the [quick start guide](http://research.libd.org/SPEAQeasy/quick-start.html).

Because SPEAQeasy is based on the [nextflow](https://www.nextflow.io/) workflow manager, it supports execution on computing clusters managed by [SLURM](https://slurm.schedmd.com/overview.html) or [SGE](https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html) without any configuration (local execution is also possible). Those with access to [docker](https://www.docker.com/) can very simply use docker containers to manage SPEAQeasy software dependencies, though we provide a script for installing dependencies for users without docker or even root privileges.


## Annotation ##

The pipeline is intended to be easily customizable regarding which annotation/ reference-related files can be used. By default, for "hg19", "hg38", and "mm10" references, the pipeline uses files provided by GENCODE; for "rn6" reference, the files are provided by Ensembl. These files are managed automatically- necessary files are downloaded whenever they are not already present, and cached for future runs.

### Configuration ###

The user can specify in the config file (determined above) the following information:

+ *Annotation version*: The variables `gencode_version_human` and `gencode_version_mouse` refer to the GENCODE release number. Similarly, the variable `ensembl_version_rat` specifies the Ensembl version for "rn6" reference.
+ *Annotation build*: The `annotation_build` variable controls whether the user wishes to include extra scaffolds. Two values are currently supported for this variable. The value "main" indicates only the main reference sequences should be included (e.g. the 25 sequences chr1-chrM for human). A value of "primary" specifies to include additional scaffolds- these definitions of "main" and "primary" come from the convention GENCODE uses in naming their reference files; however, `annotation_build` also applies to the Ensembl-based "rn6" annotation files.

### Using custom annotation ###

You may wish to provide specific reference files in place of the automatically managed files described above. In this case, you must supply the following files in the directory specified in the command-line option `--annotation [dir]`:

+ A genome assembly fasta (the reference genome to align reads to), such as the file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz). Make sure the file has the string "assembly" in the filename, to specify to the pipeline that it is the genome reference fasta.
+ Gene annotation gtf, such as the file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz)- but not gzipped. This file can have any name, so long as it ends in ".gtf".
+ A transcripts fasta, such as the file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz)- but not gzipped. Make sure to include "transcripts" anywhere in the filename (provided the file ends in ".fa") to differentiate this file from the reference genome.

*Optional files to include depending on your use-case:*

+ An ERCC index (this is a file specific to Kallisto needed for ERCC quantification, which is an optional component of the pipeline). You can find the index used by default at `[repository directory]/Annotation/ERCC/ERCC92.idx`. This file must end in ".idx".
+ A list of SNV sites at which to call variants (in .bed format). Variant calling is by default only enabled for human reference. You can find the .bed files used by default for "hg38" and "hg19" at `[repository directory]/Annotation/Genotyping/common_missense_SNVs_hg*.bed`. This file can have any name provided it has the ".bed" extension.

You must also add the `--custom_anno [label]` argument to your `run_pipeline_X.sh` script, to specify you are using custom annotation files. The "label" is a string you want to include in filenames generated from the annotation files you provided. This is intended to allow the use of potentially many different custom annotations, assigned a unique and informative name you choose each time. This can be anything except an empty string (which internally signifies not to use custom annotation).

## Manifest ##

This pipeline requires that a `samples.manifest` file exists (see the `--input` flag), to describe the samples to be processed by the pipeline. The `samples.manifest` file associates each FASTQ file with a path and ID, and allows the pipeline to automatically merge files if necessary. Each line in `samples.manifest` should have the following format:

+ *For a set of unpaired reads* `<PATH TO FASTQ FILE>(tab)<optional MD5>(tab)<sample label/id>`
+ *For paired-end sets of reads* `<PATH TO FASTQ 1>(tab)<optional MD5 1>(tab)<PATH TO FASTQ 2>(tab)<optional MD5 2>(tab)<sample label/id>`

A line of paired-end reads could look like this:

`RNA_sample1_read1.fastq    0    RNA_sample1_read2.fastq    0    sample1`

### More details regarding inputs and the manifest ###

+ The MD5(s) on each line are for compatibility with a conventional samples.manifest structure, and are not explicitly checked in the pipeline.
+ Paths must be long/full.
+ If you have a single sample split across multiple files, you can signal for the pipeline to merge these files by repeating the sample label/id on each line of files to merge.
+ Input FASTQ files can have the following file extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
+ A `samples.manifest` file cannot include both single-end and paired-end reads; separate pipeline runs should be performed for each of these read types.
+ FASTQ files must not contain "." characters before the typical extension (e.g. sample.1.fastq), since some internal functions rely on splitting file names by ".".

## Pipeline use with limited internet access ##

+ For users who do not have internet access when executing pipeline runs, you may first run `bash scripts/manual_annotation.sh`. This script must be run from the repository directory (from a machine with internet access). Modify the four lines in the "user configuration section" at the top of the script for you particular set-up. This sets up everything so that subsequent runs of the pipeline do not need an internet connection to complete.
+ Towards the end of the pipeline run, when R objects are created containing gene/exon/junction counts, some additional data is pulled from biomaRt databases by default. However, if an internet connection is not available, this extra information is not pulled (and an error occurs). Users without access to the internet during pipeline execution should provide the `--no_biomart` command flag to suppress this error and continue without the additional biomaRt data. This is the only difference between runs with and without internet access.

## Version description ##

* Version 0.8.0 (current)

    + Docker and SGE mode fully working.

    + Complete functionality for single-end type of data for human (hg19, hg38), and mouse (mm10).

        + Variant Calling 
        + Expressed Regions detection
        + Full Coverage Rdata generation
        + Transcript Counts Rdata generation


This pipeline has been successfully run in the following Operative System(s):

* [Ubuntu 16.04.4 LTS](https://www.ubuntu.com/download/alternative-downloads)


## Docker mode ##

This pipeline can run using docker containers to avoid the installation of system wide dependencies. Just follow this instructions:

* Install docker

A set of instructions for different operating systems are available on the [Docker site](https://docs.docker.com/installation/).

* Create a docker group

```bash
sudo addgroup docker
```

* Add user to docker group

```bash
sudo usermod -aG docker <your_user>
```

* Checking installation

Log out and log back in to ensure your user is running with the correct permissions.

[Test Docker installation](https://docs.docker.com/get-started/#test-docker-installation) by running:

```bash
docker run hello-world
```
You can find more information about this setup test in the [docker site](https://docs.docker.com/get-started/#test-docker-installation)

## Reference files ##

The pipeline uses many reference files during a run. Due to size limitations in git repositories, not every reference file can be versionated.

The basic Annotation and Genotyping directories are cloned with this repository. On the first run (test run or real run) for a particular species (i. e. hg38, h19, mm10 or rn6), the missing annotation files are built in a species-dependent manner.


### Genomes ###

This pipeline works for the following genomes and versions:

| Genome Build | Organism |
|--------|------|
|hg19| human |
|hg38| human |
|mm10| mouse |
|rna6| rat |


## Output data formats ##

**Variant Calling**

+ VCF files per sample, and multi-sample

**Expressed Regions detection**

+ Rdata file

**Full Coverage Rdata generation**

+ Rdata file

**Transcript Counts Rdata generation**

+ Rdata file

## Launching a real run ##

Command Example:

```
nextflow main.nf --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -profile sge
```

This command will read files from `./input`, run as __single__ end __unstranded__ samples, with __hg38__ human genome and transcriptome as the reference, __ercc__ process for spiking quantification will also run, and so will the __fullCoverage__ process to create Coverge R data. The __sge__ profile has been selected, for execution under a SGE environment.


## Authors ##

Original Pipeline

 [Emily Burke](mailto:emily.burke@libd.org>),
 [Leonardo Collado-Tores](mailto:lcolladotor@gmail.com),
 [Andrew Jaffe](mailto:andrew.jaffe@libd.org),
 [BaDoi Phan](mailto:badoi.phan@pitt.edu) 
 
Nextflow Port

 [Jacob Leonard](mailto:leonard.jacob09@gmail.com),
 [Israel Aguilar](mailto:iaguilaror@gmail.com),
 [Violeta Larios](mailto:siedracko@gmail.com),
 [Everardo Gutierrez](mailto:ever.gmillan@gmail.com)
 [Nick Eagles](mailto:nick.eagles@libd.org)

## Contact ##

* [Leonardo Collado Torres](http://lcolladotor.github.io/)
* [Winter Genomics Team](http://www.wintergenomics.com)
