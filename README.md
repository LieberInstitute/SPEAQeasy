# LIEBER INSTITUTE JAFFE LAB RNA - SEQ ANALYSIS PIPELINE #

### Summary ###

This pipeline is a RNA-seq processing tool based on the [RNAseq-pipeline](https://github.com/LieberInstitute/RNAseq-pipeline). Build on nextflow, and capable of using Docker containers and SGE task managing, this port of the RNAseq-pipeline can be used in different computer environments.

The main function of this pipeline is to produce comparable files to those used in [recount2](https://jhubiostatistics.shinyapps.io/recount/), a tool that provides gene, exon, exon-exon junction and base-pair level data.

This pipeline allows researchers to contribute data to the recount2 project even from outside the [JHPCE](https://jhpce.jhu.edu/).


### Workflow overview ###

![General Workflow](https://github.com/LieberInstitute/RNAsp/blob/master/notes/General_Workflow.png)

# Getting started #

## Software Requirements ##

+ This pipeline runs [nextflow](https://www.nextflow.io/), which requires a Java runtime. If java is not installed, you can install it on linux with `apt install default-jre`, or with a different package manager you prefer. Python 3 and pip (automatically installed with typical installations of python) are required as well.
+ Additional software configuration depends on the options available on your system/ execution environment:
    + **Using docker** (Recommended for non-JHPCE users): If docker is installed in your environment, this option requires minimal setup/ installation. From within the repository, run `bash install_software.sh "docker"`. This installs nextflow and prepares some test files, which is a one-time setup.
    + **Installing dependencies locally** (Alternative not requiring docker): The script `install_software.sh` is included in the repository, and automates the installation process. Make sure that you first just have Java (8 or later) and Python 3 installed globally (requiring root privileges). Then, from within the repository, run `bash install_software.sh "local"` for one-time setup of the pipeline.

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
|[nextflow](https://www.nextflow.io/docs/latest/getstarted.html) | >=0.27.0 (tested with 20.01.0) | `nextflow` |
|[R](https://cran.r-project.org/bin/linux/ubuntu/README.html#installation) | 3.6 | `Rscript` |
|[regtools](https://github.com/griffithlab/regtools#installation) | 0.5.1 | `regtools` |
|[RSeQC](http://rseqc.sourceforge.net/#installation) | 3.0.1 | `bam2wig.py`|
|[salmon](http://salmon.readthedocs.io/en/latest/building.html) | 1.0.0 | `salmon` |
|[samtools](http://www.htslib.org/download/) | 1.9 | `samtools` |
|[SubRead](http://bioinf.wehi.edu.au/subread-package/) | 2.0.0 | `featureCounts` |
|[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.39 | `java -jar path/to/trimmomatic.jar` |
|[wiggletools](https://github.com/Ensembl/WiggleTools) | 1.2.1 | `wiggletools` |
|[wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/) | 4 | `wigToBigWig` |

### Run on the JHPCE cluster ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. (Optional) **Adjust configuration**: hardware resource usage, software versioning, and cluster option choices are specified in *conf/jhpce.config*
3. **Modify the main script and run**: this is *run_pipeline_jhpce.sh*. The pipeline run is submitted to the cluster by executing `qsub run_pipeline_jhpce.sh`. See "Full list of command-line options" for details about modifying the script you choose.


### Run on a Sun Grid Engines (SGE) cluster ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. **Choose how to manage software dependencies**: see "Software requirements" section.
3. (Optional) **Adjust configuration**: hardware resource usage, software versioning, and cluster option choices are specified in *conf/sge.config*, if you have installed software dependencies locally, or *conf/docker_sge.config* if you will use docker.
4. **Modify the main script and run**: the main script is *run_pipeline_sge.sh*. Submit the script as a job to your cluster with `qsub run_pipeline_sge.sh`. If you are using docker, make sure to change the line `-profile sge` to `profile docker_sge`. See "Full list of command-line options" for other details about modifying the script for your use-case.

See [here](https://www.nextflow.io/docs/latest/executor.html#sge) for additional information on nextflow for SGE environments.

### Run in a SLURM environment ###

1. **Clone the repository in the current directory**: *git clone git@github.com:LieberInstitute/RNAsp.git*
2. **Choose how to manage software dependencies**: see "Software requirements" section.
3. (Optional) **Adjust configuration**: hardware resource usage, software versioning, and cluster option choices are specified in *conf/slurm.config*, if you have installed software dependencies locally, or *conf/docker_slurm.config* if you will use docker.
4. **Modify the main script and run**: the main script is *run_pipeline_slurm.sh*. Submit the script as a job to your cluster with `sbatch run_pipeline_slurm.sh`. If you are using docker, make sure to change the line `-profile slurm` to `profile docker_slurm`. See "Full list of command-line options" for other details about modifying the script for your use-case.

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
+ `--trim_mode`  Determines the conditions under which trimming occurs:
                    "skip": do not perform trimming on samples
                    "adaptive": [default] perform trimming on samples that have failed the FastQC "Adapter content" metric
                    "force": perform trimming on all samples
+ `--keep_unpaired` include this flag to keep unpaired reads output from trimming paired-end samples, for use in alignment. Default: false, as this can cause issues in downstream tools like FeatureCounts.
+ `--use_salmon`  Include this flag to quantify transcripts with Salmon rather than the default of Kallisto.
+ `--custom_anno [label]` Include this flag to indicate that the directory specified with `--annotation [dir]` includes user-provided annotation files to use instead of the default files. See the "Using custom annotation" section for more details.
+ `--force_strand` Include this flag to continue pipeline execution with a warning, when user-provided strand contrasts with inferred strandness in any sample. Default: false (halt pipeline execution with an error message if any sample appears to be a different strandness than stated by the user)
+ `--no_biomart` Include this flag to suppress potential errors in retrieving additional annotation info, such as gene symbols, from biomaRt. BiomaRt will still be queried for this info, but the pipeline will proceed without it upon failure for any reason. This option is required for users without internet access during pipeline runs.

### Nextflow Options ###

The nextflow command itself provides many additional options you may add to your "main" script. A few of the most commonly applicable ones are documented below. For a full list, type `[path to nextflow] run -h`- the full list does not appear to be documented at [nextflow's website](https://www.nextflow.io/).

+ `-w [path]` Path to the directory where nextflow will place temporary files. This directory can fill up very quickly, especially for large experiments, and so it can be useful to set this to a scratch directory or filesystem with plenty of storage capacity.
+ `-resume` Include this flag if pipeline execution halts with an error for any reason, and you wish to continue where you left off from last run. Otherwise, BY DEFAULT, NEXTFLOW WILL RESTART EXECUTION FROM THE BEGINNING.
+ `-with-report [filename]` Include this to produce an html report with execution details (such as memory usage, completion details, and much more)
+ `N [email address]` Sends email to the specified address to notify the user regarding pipeline completion. Note that nextflow relies on the `sendmail` tool for this functionality- therefore `sendmail` must be available for this option to work.

### Annotation ###

The pipeline is intended to be easily customizable regarding which annotation/ reference-related files can be used. By default, for "hg19", "hg38", and "mm10" references, the pipeline uses files provided by GENCODE; for "rn6" reference, the files are provided by Ensembl. These files are managed automatically- necessary files are downloaded whenever they are not already present, and cached for future runs.

#### Configuration ####

The user can specify in the config file (determined above) the following information:

+ *Annotation version*: The variables `gencode_version_human` and `gencode_version_mouse` refer to the GENCODE release number. Similarly, the variable `ensembl_version_rat` specifies the Ensembl version for "rn6" reference.
+ *Annotation build*: The `annotation_build` variable controls whether the user wishes to include extra scaffolds. Two values are currently supported for this variable. The value "main" indicates only the main reference sequences should be included (e.g. the 25 sequences chr1-chrM for human). A value of "primary" specifies to include additional scaffolds- these definitions of "main" and "primary" come from the convention GENCODE uses in naming their reference files; however, `annotation_build` also applies to the Ensembl-based "rn6" annotation files.

#### Using custom annotation ####

You may wish to provide specific reference files in place of the automatically managed files described above. In this case, you must supply the following files in the directory specified in the command-line option `--annotation [dir]`:

+ A genome assembly fasta (the reference genome to align reads to), such as the file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz). Make sure the file has the string "assembly" in the filename, to specify to the pipeline that it is the genome reference fasta.
+ Gene annotation gtf, such as the file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz)- but not gzipped. This file can have any name, so long as it ends in ".gtf".
+ A transcripts fasta, such as the file [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz)- but not gzipped. Make sure to include "transcripts" anywhere in the filename (provided the file ends in ".fa") to differentiate this file from the reference genome.

*Optional files to include depending on your use-case:*

+ An ERCC index (this is a file specific to Kallisto needed for ERCC quantification, which is an optional component of the pipeline). You can find the index used by default at `[repository directory]/Annotation/ERCC/ERCC92.idx`. This file must end in ".idx".
+ A list of SNV sites at which to call variants (in .bed format). Variant calling is by default only enabled for human reference. You can find the .bed files used by default for "hg38" and "hg19" at `[repository directory]/Annotation/Genotyping/common_missense_SNVs_hg*.bed`. This file can have any name provided it has the ".bed" extension.

You must also add the `--custom_anno [label]` argument to your `run_pipeline_X.sh` script, to specify you are using custom annotation files. The "label" is a string you want to include in filenames generated from the annotation files you provided. This is intended to allow the use of potentially many different custom annotations, assigned a unique and informative name you choose each time. This can be anything except an empty string (which internally signifies not to use custom annotation).

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
+ FASTQ files must not contain "." characters before the typical extension (e.g. sample.1.fastq), since some internal functions rely on splitting file names by ".".

### Pipeline use with limited internet access ###

+ For users who do not have internet access when executing pipeline runs, you may first run `bash scripts/manual_annotation.sh`. This script must be run from the repository directory (from a machine with internet access). Modify the four lines in the "user configuration section" at the top of the script for you particular set-up. This sets up everything so that subsequent runs of the pipeline do not need an internet connection to complete.
+ Towards the end of the pipeline run, when R objects are created containing gene/exon/junction counts, some additional data is pulled from biomaRt databases by default. However, if an internet connection is not available, this extra information is not pulled (and an error occurs). Users without access to the internet during pipeline execution should provide the `--no_biomart` command flag to suppress this error and continue without the additional biomaRt data. This is the only difference between runs with and without internet access.

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


#### Genomes ####

This pipeline works for the following genomes and versions:

| Genome Build | Organism |
|--------|------|
|hg19| human |
|hg38| human |
|mm10| mouse |
|rna6| rat |


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


### Run with Docker ###

* Containers

The following container versions are used in this pipeline.

| Image | Tag | Software |
|:-------------:| -----:| -----: |
| libddocker/r_3.6.1_bioc | latest | R and Bioconductor 3.10 |
| libddocker/ubuntu16.04_base | 1_v3 | Ubuntu Base |
| libddocker/infer_strandness | latest | R, Bioconductor 3.10, Kallisto |
| zlskidmore/kallisto | 0.46.0 | Kallisto |
| libddocker/fastqc | 0.11.8  | FastQC |
| libddocker/trimmomatic | 0.39 | Trimmomatic |
| libddocker/hisat2 | 2.1.0 | HISAT |
| libddocker/rseqc | 3.0.1 | RSeQC |
| libddocker/samtools | 1.9 | Samtools |
| libddocker/salmon | 0.14.1 | Salmon |
| libddocker/regtools | 0.5.1 | Regtools |
| libddocker/subread | 2.0.0 | SubRead/FeatureCounts |
| libddocker/wiggletools-1.2 | 1_v3 | Wiggletools |
| libddocker/variant_calling | 1.9 | Samtools, BCFTools |


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
