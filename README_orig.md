# LIEBER INSTITUTE JAFFE LAB RNA - SEQ ANALYSIS PIPELINE #

### Summary ###

This pipeline is a RNA-seq processing tool based on the [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline). Build on nextflow, and capable of using Docker containers and SGE task managing, this port of the RNAseq-pipeline can be used in different computer environments.

The main function of this pipeline is to produce comparable files to those used in [recount2](https://jhubiostatistics.shinyapps.io/recount/), a tool that provides gene, exon, exon-exon junction and base-pair level data.

This pipeline allows researchers to contribute data to the recount2 project even from outside the [JHPCE](https://jhpce.jhu.edu/).


### Workflow overview ###

![General Workflow](https://github.com/LieberInstitute/RNAsp/blob/master/notes/General_Workflow.png)

### Version description ###

* Version 0.8.0 (current)

    + Docker and SGE mode fully working.

    + Complete functionality for single-end type of data for human (hg19, hg38), and mouse (mm10).

        + Variant Calling 
        + Expressed Regions detection
        + Full Coverage Rdata generation
        + Transcript Counts Rdata generation


### Installation ###

##### Working OS #####

This pipeline has been successfully run in the following Operative System(s):

* [Ubuntu 16.04.4 LTS](https://www.ubuntu.com/download/alternative-downloads)

**Note on dependencies:** This pipeline can run in System mode (using tools installed in the system), 
or in Docker mode (using docker to handle software dependencies). Install dependencies accordingly.

##### System Mode #####

If you are going to run this pipeline without Docker mode enabled, please verify that your system has the following tools and versions:

Software | Version | Command used by the pipeline |
|:-------------:| -----:| -----: |
|[bcftools](http://www.htslib.org/download/) | 1.6 | `bcftools` |
|[fastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc) | 0.11.4 | `fastqc` |
|[hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml#obtaining-hisat2) | 2.0.4 | `hisat2`, `hisat2-build` |
|[htslib](http://www.htslib.org/download/) | 1.2 | `tabix` |
|[java](http://www.oracle.com/technetwork/java/javase/downloads/index.html) | 8 | `java` |
|[kallisto](https://pachterlab.github.io/kallisto/source) | 0.44.0 | `kallisto` |
|[nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | 0.29.1.4804 | `nextflow` |
|[R](https://cran.r-project.org/bin/linux/ubuntu/README.html#installation) | 3.4.4 | `Rscript` |
|[regtools](https://github.com/griffithlab/regtools#installation) | 0.5.0 | `regtools` |
|[RSeQC](http://rseqc.sourceforge.net/#installation) | 2.6.4 | `infer_experiment.py`, `bam2wig.py`|
|[salmon](http://salmon.readthedocs.io/en/latest/building.html) | 0.9.1 | `salmon` |
|[samtools](http://www.htslib.org/download/) | 1.2 | `samtools` |
|[SubRead](http://bioinf.wehi.edu.au/subread-package/) | 1.6.0 | `featureCounts` |
|[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.36 | `java -jar path/to/trimmomatic.jar` |
|[vcf-tools](https://vcftools.github.io/index.html) | 0.1.15 | `vcf-merge` |
|[wiggletools](https://github.com/Ensembl/WiggleTools) | 1.2 | `wiggletools` |
|[wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | 4 | `wigToBigWig` |

##### Docker mode #####

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

##### SGE mode #####

This pipeline can be scaled up using a Sun Grid Engine cluster or a compatible platform (Open Grid Engine, Univa Grid Engine, etc), as described [here](https://www.nextflow.io/docs/latest/executor.html#sge).

Contact your computing cluster administration to ensure you have access to a queue for submiting jobs, and for information about the best practices for resource management during pipeline execution.

### Pipeline setup ###

Once the previous dependencies have been met, clone this repository via:

````
        git clone https://github.com/LieberInstitute/RNAsp.git
````

Continue with the Configuration step in preparation for a test run.

### Configuration ###

Thanks to the nextflow framework this pipeline can run on your local machine, or in a SGE cluster; with or without using docker containers.

But first, you need to configure some variables in the following files:

* **conf/command_paths.config**: this file defines the paths used by the pipeline to make some required command calls.  

    + This conf file can allow you to test the pipeline even if some dependencies are not globally installed or available on the PATH.
    + **Important**: change the values to match your system environment.

* **conf/mem.config**: this file defines the amount of computational resources that Nextflow requests for every process. This conf file is used by default if no other `-profile` is requested.

    + By default, this configuration file assumes the local environment has 16G of memory and 8 CPUs

* **conf/sge.config**: this file defines variables used by SGE during job submitions, mainly computational resources.

    + **Important**: change the ***queue*** variable to a valid queue where your user is allowed to submit jobs.
    + **Important**: change the ***penv*** variable to a valid parallel environment according to your cluster setup.
    + **Important**: this conf file is not used by default. It must be requested by using the `-profile sge` option of the nextflow command.

* **conf/sge_large.config**: see **sge.config**, this file is similar but should be configured to request heavier use of computational resources.

    + **Important**: this conf file is not used by default. It must be requested by using the `-profile sge_large` option of the nextflow command.
    
* **conf/jhpce.config**: this file is designed for use at the JHPCE cluster, the SGE cluster on which some of the pipeline was developed. It specifies computational resources and cluster options intended for this environment; use sge.config or sge_large.config for settings that apply more generally to SGE clusters.

### Test run ###

* System mode simple test

After proper configuration has been made in the _**conf/command_paths.config**_ file, you can launch a test by executing:

````
bash run_test_system.sh
````

This will launch a **local run** of the complete pipeline.

* System + SGE simple test

After proper configuration has been made in the _**conf/command_paths.config**_ AND the _**conf/sge.config**_ files, launch this test by executing:

````
bash run_test_sge.sh
````

* Docker simple test

**NOTE**: First read the _**Run with Docker**_ section of this readme.

Launch this test by executing:

````
bash run_test_docker.sh
````

### Reference files ###

The pipeline uses many reference files during a run. Due to size limitations in git repositories, not every reference file can be versionated.

The basic Annotation and Genotyping directories are cloned with this repository. On the first run (test run or real run) for a particular species (i. e. hg38, h19, mm10 or rn6), the missing annotation files are built in a species-dependent manner.

A tree view for full Annotation and Genotyping directories can be consulted in ***notes/reference_directories_structure.md***.

#### Genomes ####

This pipeline works for the following genomes and versions:

| Genome Build | Organism |
|--------|------|
|hg19| human |
|hg38| human |
|mm10| mouse |
|rna6| rat |

### Input data formats ###

This RNA-Seq pipeline can handle several combinations of single/paired read, human/mouse/rat reference, and unstranded'forward/reverse options. The following __mandatory__ options are available:

__reference__

````
	Human:
	--reference "hg38"
	  OR
	--reference "hg19"

	Mouse:
	--reference "mm10"

	Rat:
	--reference "rn6"
````

__sample__

````
	--sample "single"
	  OR
	--sample "paired"
````

__stranded__

````
	Stranded:
	--strand "forward"
	  OR
	--strand "reverse"

	Unstranded:
	--strand "unstranded"
````


### File Naming ###

NOTICE: File names must not contain "." in the name because the pipeline operates on file names by splitting along the "." to determine prefixes. If needed, change all "." characters to "\_"

### Merging Required

The pipeline can handle sample fastq.gz files that require merging. Files that need merging should be named in the following format:


__single__

````
	Sample part 1: "{prefix}_read1.fastq.gz"
	Sample part 2: "{prefix}_read2.fastq.gz"
````

__paired__

````
	First in pair, Sample part 1: "{prefix}_1_read1.fastq.gz"
	First in pair, Sample part 2: "{prefix}_1_read2.fastq.gz"
	Second in pair, Sample part 1: "{prefix}_2_read1.fastq.gz"
	Second in pair, Sample part 2: "{prefix}_2_read2.fastq.gz"
````

### Output data formats ###

**Variant Calling**

+ VCF files per sample, and multi-sample

**Expressed Regions detection**

+ Rdata file

**Full Coverage Rdata generation**

+ Rdata file

**Transcript Counts Rdata generation**

+ Rdata file

### Launching a real run ###

Command Example:

```
nextflow main.nf --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -profile sge
```

This command will read files from `./input`, run as __single__ end __unstranded__ samples, with __hg38__ human genome and transcriptome as the reference, __ercc__ process for spiking quantification will also run, and so will the __fullCoverage__ process to create Coverge R data. The __sge__ profile has been selected, for execution under a SGE environment.

### Email notifications

We use the built-in notification system form nextflow, as described [here](https://www.nextflow.io/docs/latest/mail.html?highlight=email#workflow-notification) :

> Nextflow includes a built-in workflow notification features that automatically sends a notification message when a workflow execution terminates.   
To enable simply specify the -N option when launching the pipeline execution. For example:  

````
nextflow run main.nf <pipeline options> -N <recipient email address>
````  

This will send a notification mail when the execution completes.  

**Warning**: By default the notification message is sent by using the `sendmail` system tool which is assumed to be available in the computer where Nextflow is running. Make sure it's properly installed and configured.

### Run with Docker ###

* Containers

The following container versions are used in this pipeline.

| Container | tag | software |
|:-------------:| -----:| -----: |
| r3.4.3_base | 1_v3 | R Base |
| ubuntu16.04_base | 1_v3 | Ubuntu Base |
| kallisto_v0.43.1 | 1_v3 | Kallisto |
| fastqc_v0.11.5 | 1_v3 | FastQC |
| trimmomatic-0.36 | 1_v3 | Trimmomatic |
| hisat2-2.0.4 | 1_v3 | HISAT |
| rseqc-v2.6.4 | 1_v3 | RSeQC |
| samtools-1.3.1 | 1_v3 | Samtools |
| salmon-0.9.1 | 1_v3 | Salmon |
| regtools-0.3.0 | 1_v3 | Regtools |
| subread-1.6.0 | 1_v3 | SubRead |
| wiggletools-1.2 | 1_v3 | Wiggletools |
| bcftools-1.3.1 | 1_v3 | BCFTools |
| vcftools-v0.1.15 | 1_v3 | VCFTOols |

Once Docker is installed, required docker images can be pulled down from **[DockerHub](https://hub.docker.com/u/libdocker/)**, 
or built locally. The dockerfiles for each container can be found in the ./dockerfiles folder, and are built 
using versioning to maintain reproducibility, based on [best practices for writting dockers](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/).

* Docker Images

Required docker containers can be obtained locally in 1 of 2 ways:

1) Docker pull command

```
docker pull libdocker/${CONTAINER}
```

2) Docker build command

```
docker build -t libdocker/${CONTAINER} .
```	

An example of how to build and deploy all the required docker images can be found in _dockerfiles/make.log_

* Potential Issues with Docker

R-Base Dockerfile: This dockerfile holds all of the R packages needed throughout the pipeline. 
In order to maintain a specific version combination for required software, all packages are manually installed 
from source. However, this causes the dockerfile to have numerous commit layers, and risks reaching maximum depth ~125 
layers. In this case, packages can be consolidated and installed by name (rather than source link) to decrease the 
number of layers. OR, the existing docker image can be squashed using --squash (with docker in --experimental mode) to 
flatten the existing image

The dockerfiles declare specific versions of the required software, since software becomes outdated and links need to 
be updated, the dockerfiles will become less likely to build successfully. In this case its much easier to simply run 
make pull in the dockerfiles directory to pull the current working version of each image.

### Run with SGE ###

For scalability, this pipeline uses the executor component from Nextflow, as described [here](https://www.nextflow.io/docs/latest/executor.html); especifically, we use the [SGE](https://www.nextflow.io/docs/latest/executor.html#sge) integration capabilities to manage process distribution and computational resources.

The _conf/sge.config_ and _conf/sge_large.config_ must be properly configured before launching SGE runs. Said configuration files define variables regarding queue, parallelization environments and resources requested by every process in the pipeline. This allows the fine tunning of resource consumption.

### Authors ###

Original Pipeline

 [Emily Burke](mailto:emily.burke@libd.org>),
 [Leonardo Collado-Tores](mailto:lcolladotor@gmail.com),
 [Andrew Jaffee](mailto:andrew.jaffe@libd.org),
 [BaDoi Phan](mailto:badoi.phan@pitt.edu) 
 
Nextflow Port

 [Jacob Leonard](mailto:leonard.jacob09@gmail.com),
 [Israel Aguilar](mailto:iaguilaror@gmail.com),
 [Violeta Larios](mailto:siedracko@gmail.com),
 [Everardo Gutierrez](mailto:ever.gmillan@gmail.com)

### Contact ###

* [Leonardo Collado Torres](http://lcolladotor.github.io/)
* [Winter Genomics Team](http://www.wintergenomics.com)