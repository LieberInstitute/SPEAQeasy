# LIEBER INSTITUTE JAFFE LAB RNA - SEQ ANALYSIS PIPELINE #

### Summary ###

Aenean sodales velit at elementum blandit. Donec lobortis tempor aliquam. Maecenas tempor egestas ipsum eget congue. In arcu magna, venenatis efficitur sapien nec, feugiat mattis urna...

_##**TODO**##_: Brief description of the pipeline

_##**TODO**##_: Show a simplified workflow of the pipeline, from a notes/*.png image (done in draw.io or something similar)

### Version description ###

* Version 0.7.5 (current)

    + **Test run validation**:  for hg38, hg10, mm10, and rn6; with the command:  
`--small_test --sample "single" --strand "unstranded" --ercc -with-report -with-dag -N user@email.com`
    + **Test run validation**: --fullCov option working for hg38, hg19 and mm10; _rn6 **requires debugging in create_count_objects-rat.R** script_.
    + **Process validation**: All procceses for Annotation references construction, validated for hg38, hg19, mm10, and rn6.
    + **Portability feature**: added conf/command paths.config file for defining paths to commands and essential .py scripts.
    + **Basic feature**: --ercc and --fullcov options functional.
    + **Documentation expansion**: reestructured README.md; added basic dependencies info; added test run instructions; added email notification info; added Reference files info; added notes on Reference file directories; added some configuration info; added author info.

### Installation ###

##### Working OS #####

This pipeline has been successfully run in the following Operative System(s):

* [Ubuntu 16.04.4 LTS](https://www.ubuntu.com/download/alternative-downloads)

##### Dependencies: Bioinformatics software #####

Please verify that your system has the following tools and versions:

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

_##**TODO**##_: Add optional dependencies for SGE


_##**TODO**##_: Add optional dependencies for docker

Install docker

A set of instructions for different operating systems are available on the Docker link **(https://docs.docker.com/installation/).

**Create a docker group**

```bash
sudo addgroup docker
```

**Add user to docker group**  

Here user is ``ec2-user``

```bash
sudo usermod -aG docker ec2-user
```
Log out and log back in.

This ensures your user is running with the correct permissions.

Verify your work by running ``docker`` without ``sudo``.

```bash
docker run hello-world
```
## Docker Security...

This post reviews the various security implications of using Docker to run applications within containers, and how to address them:
[How Secure are Containers?](https://blog.docker.com/2013/08/containers-docker-how-secure-are-they/#more-697)

> Docker containers are, by default, quite secure; especially if you take care of running your processes inside the containers as non-privileged users (i.e. non root).


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
    + **Important**: change the values to match your system environment.
    + This conf file can allow you to test the pipeline even if some dependencies are not globally installed or available on the PATH.

* **conf/sge.config**: this file defines variables used by SGE during job submitions.
    + **Important**: change the ***queue*** variable to a valid queue where your user is allowed to submit jobs.
    + **Important**: change the ***penv*** variable to a valid parallel environment according to your cluster setup.

.

_##**TODO**##_: describe conf/docker.config

_##**TODO**##_: describe conf/mem.config

_##**TODO**##_: describe conf/sge_large.config

### Test run ###

* Local simple test

After proper configuration has been made in the _**conf/command_paths.config**_ file, you can launch a test by executing:

````
bash run_test_simple.sh [-N user@email.com]
````

This will launch a **local run** of the complete pipeline.

* SGE simple test

After proper configuration has been made in the _**conf/command_paths.config**_ AND the _**conf/sge.config**_ files, launch this test by executing:

````
bash run_test_sge.sh [-N user@email.com]
````

* Docker simple test (**NOT IMPLEMENTED YET**)

**NOTE**: First read the _**Run with Docker**_ section of this readme.


After propper configuration has been made in the _**conf/command_paths.config**_ AND the _**conf/docker.config**_ files, launch this test by executing:

````
bash run_test_docker.sh [-N user@email.com]
````

### Reference files ###

The pipeline uses many reference files during a run. Due to size limitations in git repositories, not every reference file can be versionated.

The basic Annotation and Genotyping directories are cloned with this repository. On the first run (test run or real run) for a particular species (i. e. hg38, h19, mm10 or rn6), the missing annotation files are built in a species-dependent manner.

A tree view for full Annotation and Genotyping directories can be consulted in ***notes/reference_directories_structure.md***.


### Process description ###

Sed dictum tristique bibendum. Nulla posuere lacus nec auctor consequat. Ut a sodales orci.
 
_##**TODO**##_: Describe for every process step of the pipeline, what it does (order it accordingly to what it is described in the main.nf header)

### Input data formats ###

Sed bibendum felis eu consequat aliquet. Donec elementum rhoncus massa, et egestas tortor condimentum volutpat. Nam nunc sapien, laoreet quis pulvinar in, finibus vel mauris. Etiam et tellus ligula...

_##**TODO**##_: describe species and data types (single, paired, etc.) accepted by the pipeline.

## Genomes

| Genome Build | Organism |
|--------|------|
|hg19| human |
|hg38| human |
|mm10| mouse |
|rna6| rat |

_##**TODO**##_: Make notes about file naming, for normal runs, and for --merged runs

_##**TODO**##_: Make notes about disabled modules for mm10 and rn6, if any

### Output data formats ###

Quisque vitae venenatis lorem. Nulla id dui euismod, semper ipsum a, auctor magna. Fusce eget feugiat augue, ut mattis felis...

_##**TODO**##_: describe the many output files produced by the pipeline.

_##**TODO**##_: Consult with Lieber which files are final outputs and which are temporary files

### Launching a real run ###

Suspendisse porttitor, nibh id euismod consectetur, lectus nisl posuere nisi, et egestas dui tellus quis dui. Suspendisse dignissim justo ac aliquam efficitur...

_##**TODO**##_: describe how to launch a normal run

_##**TODO**##_: describe the many options of the pipeline, the flags and what they mean

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

Duis imperdiet, nisl ac imperdiet vehicula, ipsum arcu dictum justo, vitae sollicitudin tellus tortor a ligula. Curabitur id sapien faucibus, luctus neque vitae, convallis purus...

_##**TODO**##_: describe how docker integration works in this pipeline

_##**TODO**##_: describe how to build or pull the dockers.

Once Docker is installed,images can be pulled down from **[DockerHub] (https://hub.docker.com/u/libdocker/)**. The dockerfiles for each container can be found in the ./dockerfiles folder, and are built using versioning to maintain reproducibility and base on best practices in writting dockers ** (https://docs.docker.com/develop/develop-images/dockerfile_best-practices/)**. Docker Images can be obtained locally in 1 of 2 ways.

* Containers

The following conatainers versions are used in this pipeline. For the Ubuntu Base and R Base images, look to the dockerfiles to see specifics on installed packages.

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

Docker pull command

```
docker pull libdocker/${CONTAINER}
```

Docker build command

```
docker build -t libdocker/${CONTAINER} .
```	

An example build and deploy of all the images can be found in ./dockerfiles/make.log

Potential Issues

R-Base Dockerfile: This dockerfile holds all of the R packages needed throughout the pipeline. In order to maintain a specific version combination for required software, all packages are manually installed from source. However, this causes the dockerfile to have numerous commit layers, and risks reaching maximum depth ~125 layers. In this case, packages can be consolidated and installed by name (rather than source link) to decrease the number of layers. OR, the existing docker image can be squashed using --squash (with docker in --experimental mode) to flatten the existing image

The dockerfiles declare specific versions of the required software, and as software becomes outdated and links need to be updated, the dockerfiles will become less likely to build successfully. In this case its much easier to simply run make pull in the dockerfiles directory to pull the current working version of each image.

_##**TODO**##_: describe the --with-docker flag

### Run with SGE ###

Nulla ultrices ligula et nunc laoreet pretium. Suspendisse placerat sapien velit, a vulputate justo volutpat sit amet. Praesent massa dui, varius id sodales a, maximus eget tortor...

_##**TODO**##_: describe how SGE integration works in this pipeline

_##**TODO**##_: describe how to configure the cong/sge* configuration files, regarding queue, profile and resources request

### Pipeline directory structure ###

Proin euismod ligula ac est sagittis, ac egestas ante malesuada. Nam hendrerit dui eu nunc molestie maximus. Aliquam faucibus sapien eget ante cursus, non accumsan ligula tincidunt...

_##**TODO**##_: make a tree view of a final directory print from all the sucessfull runs.

_##**TODO**##_: describe in brief wvery file in the tree

### Authors ###

_##**TODO**##_: complete email adresses for the Lieber team

Original Pipeline

 [Emily Burke](mailto:user@email.com>),
 [Leonardo Collado-Tores](mailto:fellgernon@gmail.com),
 [Andrew Jaffee](mailto:user@email.com),
 [BaDoi Phan](mailto:user@email.com) 
 
Nextflow Port

 [Jacob Leonard](mailto:leonard.jacob09@gmail.com),
 [Israel Aguilar](mailto:iaguilaror@gmail.com),
 [Violeta Larios](mailto:siedracko@gmail.com),
 [Everardo Gutierrez](mailto:ever.gmillan@gmail.com)

### Contact ###

* [Leonardo Collado Torres](http://lcolladotor.github.io/)
* [Winter Genomics Team](http://www.wintergenomics.com)
