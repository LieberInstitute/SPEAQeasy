# LIEBER INSTITUTE JAFFE LAB RNA - SEQ ANALYSIS PIPELINE #

### Summary ###

This pipeline is a RNA-seq processing tool based on the [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline). Build on nextflow, and capable of using Docker containers and SGE task managing, this port of the RNAseq-pipeline can be used in different computer environments.

The main function of this pipeline is to produce comparable files to those used in [recount2](https://jhubiostatistics.shinyapps.io/recount/), a tool that provides gene, exon, exon-exon junction and base-pair level data.

This pipeline allows researchers to contribute data to the recount2 project even from outside the [JHPCE](https://jhpce.jhu.edu/).


### Workflow overview ###

![General Workflow](https://github.com/LieberInstitute/RNAsp/blob/master/notes/General_Workflow.png)

# Getting started #

## Software Requirements ##

+ This pipeline runs [nextflow](https://www.nextflow.io/), which requires a Java runtime. If java is not installed, you can install it on linux with `apt install default-jre`, or with a different package manager you prefer. Python 2.7 is required as well.
+ Additional software configuration depends on the options available on your system/ execution environment:
    + **Using docker** (Recommended for non-JHPCE users): If docker is installed in your environment, this option requires minimal setup/ installation. From within the repository, run `bash install_software.sh "docker"`. This installs nextflow and prepares some test files, which is a one-time setup.
    + **Installing dependencies locally** (Alternative not requiring docker): The script `install_software.sh` is included in the repository, and automates the installation process. Make sure that you first just have Java (8 or later) and Python 2.7 installed globally (requiring root privileges). Then, from within the repository, run `bash install_software.sh "local"` for one-time setup of the pipeline.

### Advanced info regarding installation ###

+ The script `install_software.sh` need only be run once. If you are installing software to run the pipeline locally, all dependencies are installed into `[repo directory]/Software/`, and `[repo directory]/conf/command_paths_long.config` is configured to show nextflow the default installation locations of each software tool. Thus, this config file can be tweaked to manually point to different paths, if need be (though this shouldn't be necessary).
+ Nextflow supports the use of Lmod modules to conveniently point the pipeline to the bioinformatics software it needs. If you neither wish to use docker nor wish to install the many dependencies locally-- and already have Lmod modules on your cluster-- this is another option. In the appropriate config file (as determined in step 3 in the section you choose below), you can include a module specification line in the associated process (such as `module = 'hisat2/2.1.0'` for buildHISATindex) as configured in *conf/jhpce.config*. In most cases this will be more work to fully configure, and so running the pipeline with docker or locally installing software is generally recommended instead. See [nextflow modules](https://www.nextflow.io/docs/latest/process.html#module) for some more information.

Here is the full list of software used by this pipeline:
    
Software | Version | Command used by the pipeline |
|:-------------:| -----:| -----: |
|[bcftools](http://www.htslib.org/download/) | 1.9 | `bcftools` |
|[fastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc) | 0.11.8 | `fastqc` |
|[hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml#obtaining-hisat2) | 2.1.0 | `hisat2`, `hisat2-build` |
|[htslib](http://www.htslib.org/download/) | 1.9 | `tabix` |
|[java](http://www.oracle.com/technetwork/java/javase/downloads/index.html) | 8 | `java` |
|[kallisto](https://pachterlab.github.io/kallisto/source) | 0.46.1 | `kallisto` |
|[nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | >=0.27.0 (tested with 18.10.0) | `nextflow` |
|[R](https://cran.r-project.org/bin/linux/ubuntu/README.html#installation) | 3.6 | `Rscript` |
|[regtools](https://github.com/griffithlab/regtools#installation) | 0.5.1 | `regtools` |
|[RSeQC](http://rseqc.sourceforge.net/#installation) | 2.6.4 | `bam2wig.py`|
|[salmon](http://salmon.readthedocs.io/en/latest/building.html) | 1.0.0 | `salmon` |
|[samtools](http://www.htslib.org/download/) | 1.9 | `samtools` |
|[SubRead](http://bioinf.wehi.edu.au/subread-package/) | 2.0.0 | `featureCounts` |
|[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.39 | `java -jar path/to/trimmomatic.jar` |
|[wiggletools](https://github.com/Ensembl/WiggleTools) | 1.2.1 | `wiggletools` |
|[wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | 4 | `wigToBigWig` |

### Run on the JHPCE cluster ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. (Optional) **Adjust configuration**: hardware resource usage, software versioning, and cluster option choices are specified in *conf/jhpce.config*
3. **Modify the main script and run**: this is *run_pipeline_jhpce_qsub.sh*. The pipeline run is submitted to the cluster by executing `qsub run_pipeline_jhpce_qsub.sh`. Alternatively, you may run the pipeline interactively via `bash run_pipeline_jhpce.sh`. See "Full list of command-line options" for details about modifying the script you choose.


### Run on a Sun Grid Engines (SGE) cluster ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. **Choose how to manage software dependencies**: see "Software requirements" section.
3. (Optional) **Adjust configuration**: hardware resource usage, software versioning, and cluster option choices are specified in *conf/sge.config*, if you have installed software dependencies locally, or *conf/docker_sge.config* if you will use docker.
4. **Modify the main script and run**: the main script is *run_pipeline_sge.sh*. Run the pipeline interactively with `bash run_pipeline_sge.sh`, or submit as a job to your cluster with `qsub run_pipeline_sge.sh`. If you are using docker, make sure to change the line `-profile sge` to `profile docker_sge`. See "Full list of command-line options" for other details about modifying the script for your use-case.

See [here](https://www.nextflow.io/docs/latest/executor.html#sge) for additional information on nextflow for SGE environments.

### Run in a SLURM environment ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. **Choose how to manage software dependencies**: see "Software requirements" section.
3. (Optional) **Adjust configuration**: hardware resource usage, software versioning, and cluster option choices are specified in *conf/slurm.config*, if you have installed software dependencies locally, or *conf/docker_slurm.config* if you will use docker.
4. **Modify the main script and run**: the main script is *run_pipeline_slurm.sh*. Run the pipeline interactively with `bash run_pipeline_slurm.sh`, or submit as a job to your cluster with `sbatch run_pipeline_slurm.sh`. If you are using docker, make sure to change the line `-profile slurm` to `profile docker_slurm`. See "Full list of command-line options" for other details about modifying the script for your use-case.

### Run locally ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. **Choose how to manage software dependencies**: see "Software requirements" section.
3. (Optional) **Adjust configuration**: hardware resource usage and other configurables are located in *conf/local.config*, if you have installed software dependencies locally, or *conf/docker_local.config* if you will use docker. Note that defaults assume access to 8 CPUs and 16GB of RAM.
4. **Modify the main script and run**: the main script is *run_pipeline_local.sh*. If you are using docker, make sure to change the line `-profile local` to `profile docker_local`. After configuring options for your use-case (See "Full list of command-line options"), simply run on the command-line with `bash run_pipeline_local.sh`.
  
Note that the configuration files also include command-line options passed to many of the software tools (such as minimum mapping quality used in samtools for filtering). This gives control over many of the parameters in the pipeline that we deemed to involve preference, or to involve variability among use-cases.

## Sharing the pipeline among many users ##
  
+ A common use-case may involve wanting to set up this pipeline once, and have potentially many users running this pipeline without additional work from the users. This can be achieved by following the above procedure to first set up the pipeline; any user wishing to execute the pipeline from a different location may copy the "main" script (`run_pipeline_[executor].sh`). Then, in this copy of the "main" script, simply change `nextflow main.nf` to `[path to nextflow executable] [path to main.nf in the original repository]` to complete the setup.

## Full list of command-line options ##

### Mandatory Parameters ###

+ `--sample`		"single" or "paired": the orientation of your reads
+ `--strand`		"unstranded", "forward, or "reverse": the strandness of your reads. Since strandness is inferred by sample in the pipeline, this option informs the pipeline to generate appropriate warnings if unexpected strandness is inferred.
+ `--reference`		"hg38", "hg19", "mm10", or "rn6": the reference genome to which reads are aligned

### Optional Parameters ###

+ `--experiment`	Name of the experiment being run (ex: "alzheimer"). Defaults to "Jlab_experiment"
+ `--prefix`	Defines the prefix of the input files (not used to detect files)
+ `--input`		The path to the directory with the "samples.manifest" file. Defaults to "./input" (relative to the repository)
+ `--output`  The path to the directory to store pipeline output files/ objects. Defaults to "./results" (relative to the repository)
+ `--unalign`		Include this flag to save discordant reads after the alignment step (false/ not included by default)
+ `--annotation`	The path to the directory containing pipeline annotations. Defaults to "./Annotations" (relative to the repository). If annotations are not found here, the pipeline includes a step to build them.
+ `--ercc`			Include this flag to enable ERCC quantification with Kallisto
+ `--fullCov`		Flag to perform full coverage analysis
+ `--small_test`	Uses sample files located in the test folder as input. Overrides the "--input" option.
+ `--force_trim`  Include this flag to perform triming on all inputs. By default, only inputs failing fastQC on the adapter content metric are trimmed.
+ `--use_salmon`  Include this flag to quantify transcripts with Salmon rather than the default of Kallisto.

### Manifest ###

This pipeline requires that a `samples.manifest` file exists (see the `--input` flag), to describe the samples to be processed by the pipeline. The `samples.manifest` file associates each FASTQ file with a path and ID, and allows the pipeline to automatically merge files if necessary. Each line in `samples.manifest` should have the following format:

+ *For a set of unpaired reads* `<PATH TO FASTQ FILE>(tab)<optional MD5>(tab)<sample label/id>`
+ *For paired-end sets of reads* `<PATH TO FASTQ 1>(tab)<optional MD5 1>(tab)<PATH TO FASTQ 2>(tab)<optional MD5 2>(tab)<sample label/id>`

A line of paired-end reads could look like this:

`RNA_sample1_read1.fastq    0    RNA_sample1_read2.fastq    0    sample1`

#### More details regarding inputs and the manifest ####

+ The MD5(s) on each line are for compatibility with a conventional samples.manifest structure, and are not explicitly checked in the pipeline.
+ Paths must be long/full.
+ If you have a single sample split across multiple files, you can signal for the pipeline to merge these files by repeating the sample label/id on each line of files to merge.
+ Input FASTQ files can have the following file extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
+ A `samples.manifest` file cannot include both single-end and paired-end reads; separate pipeline runs should be performed for each of these read types.


### Version description ###

* Version 0.8.0 (current)

    + Docker and SGE mode fully working.

    + Complete functionality for single-end type of data for human (hg19, hg38), and mouse (mm10).

        + Variant Calling 
        + Expressed Regions detection
        + Full Coverage Rdata generation
        + Transcript Counts Rdata generation


This pipeline has been successfully run in the following Operative System(s):

* [Ubuntu 16.04.4 LTS](https://www.ubuntu.com/download/alternative-downloads)


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


### File Naming ###

NOTICE: File names must not contain "." in the name because the pipeline operates on file names by splitting along the "." to determine prefixes. If needed, change all "." characters to "\_"


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


### Authors ###

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

### Contact ###

* [Leonardo Collado Torres](http://lcolladotor.github.io/)
* [Winter Genomics Team](http://www.wintergenomics.com)
