MACHUAHUITL RNA-Seq Analysis Pipeline
====================================

## Nextflow 

This RNA-Seq pipeline runs on the distributed processing framework, [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html). 

	Installation Instructions:

	```
	Nextflow is distributed as a self-contained executable package, which means that it does not
	require any special installation procedure.

	It only needs two easy steps:

	  1. Download the executable package by copying and pasting the following command in your
	  terminal window: wget -qO- https://get.nextflow.io | bash. It will create the nextflow
	  main executable file in the current directory.
	  2. Optionally, move the nextflow file to a directory accessible by your $PATH variable
	  (this is only required to avoid remembering and typing the full path to nextflow each time
	  you need to run it).
	```

The pipeline works by reading a main.nf file, and using the nexflow.config file to determine the execution and processing environment. Additionally, the nextflow.config file decalres logging parameters, and points to the execution specific config files located in the '''./conf''' folder. The development profile is set to ```standard``` and uses __docker__ to execute and maintain the required software versions. __BEFORE EXECUTING: The development enviornment and subsequent memory and cpu requirements defined for each task is held in the mem.config file, and is declared for an environment with the available resources: 8 Cores, 16 GB Memory__ The pipeline works by declaring parameters in the form ```params.parameter_name``` that can be used as command line options. These parameters can be used to define the pathway of the pipeline, and declare subsequent inputs and values. Files in this pipeline are ingested through input files and values in the form of __channels and declarations__.


__Process Overview__

	```
	Channels/Declarations(1) -> Process(1) -> Channel(2)
	Channel(2) + Channels/Declarations(3)-> Process(2)
	```
	Where channels represent the flow of input and output files in a list ingestion structure.


__Channels__

	```
	This method operates by generating lists of values in the form of integers, values and
	files. Channels can be manipulated by operators to form single value channels, paired
	value channels, prefix grouped channels, and many more forms of aggregation into the
	processing pipeline. These values can range from integers and strings, to files and lists of
	files grouped by similar name. Nextflow creates a reverse dependecy map for the whole
	pipeline to determine where to expect file outputs, where outputs from previous tasks need
	to be ingested into consecutive processes, as well as where to look for existing files for
	ingestion. Channels can be manipulated with operators that are used to group/split/
	merge/etc. the channels into different sized lists and values. This is the best way to
	ingest multiple files into a single process for scenarios such as sample aggregations and
	index ingestion. For examples, look to the pipeline_breakdown.txt for more details on how
	channels, operators work and processes.
	```

__Declarations__

	```
	This method operates simply by declaring a variable to be equal to a value such as a string
	of integer, or declaring it as a file() with a path. This method acts as a constant
	input, and can be reused by multiple processes without re-declaration. Look to the
	main.nf file starting at line 437 to see reference specific file inputs and parameter
	declarations.
	```


### Process

Processes are the content and commands of the pipeline, and represent the actual execution of pipeline steps. After a process is declared, the following segments define the way information is handled:

__echo__

	```
	Boolean used to echo the command processes to the log
	Ex:
	    echo true
	```

__tag__

	```
	String passed to the process to be used as the process label when a process is executed
	Ex: 
	    tag "Prefix: $untrimmed_prefix | Sample: [ $fastqc_untrimmed_input ]"
	```

__publishDir__

	```
	The publishDir is used to determine where to copy or move the output files from the process
	Ex: 
	    publishDir "${params.basedir}/FastQC/Untrimmed",mode:'copy'
	```

__input__

	```
	Statements that ingest files from channels or declarations
	Ex:
	    set val(samples_prefix), file(manifest_samples) from manifest_creation
	    file hisat_index from hisat_index
	```

__output__

	```
	Statements that define the output files into channels, and will be copied to the publish directory
	Ex:
	    file "*"
	    set val("${sam_to_bam_prefix}"), file("${sam_to_bam_prefix}*.sorted.bam"), file("${sam_to_bam_prefix}*.sorted.bam.bai") into infer_experiment_inputs
	```

__script__

	```
	Using script to define the process environment creates a sudo-bash script that is executed in the following way
	Ex:
	    script:
	    original_bam = "${sam_to_bam_prefix}_accepted_hits.bam"
	    sorted_bam = "${sam_to_bam_prefix}_accepted_hits.sorted"
	    samtobam_cores = "${params.samtobam_cores}"
	    """
	    samtools view -bh -F 4 $sam_to_bam_input > $original_bam
	    samtools sort -@ $samtobam_cores $original_bam -o ${sorted_bam}.bam
	    samtools index ${sorted_bam}.bam
	    """
	```

__shell__

	```
	Using shell to define the process environment creates a sudo-shell to execute commands and custom scripts
	Ex:
	    shell:
	    '''
	    python /usr/local/bin/infer_experiment.py \
	    -i !{bam_file} \
	    -r !{bedfile} \
	    1> !{infer_prefix}_experiment.txt \
	    2> !{infer_prefix}_experiment_summary_out.txt
	    '''
	```

__workdir__

	```
	The work directory is the place in which Nextflow saves the input files and values, runs the
	individual process on the sample files, and then captures the output files into channels and
	the publishDir. This directory will take up a lot of space as it copies all files needed
	each time it is run. The workDir can be changed with a flag on the nextflow command:
	Ex:
	    -w "/dump/work"
	```


## Docker

The development version of this pipeline is written using docker as the local execution environmnet. Install Docker in your local environmnet if you'd like to use this feature.

	__Version Info__ (during development)

	```
	Client:
	  Version:	17.12.1-ce
	    API version:	1.35
	    Go version:	go1.9.4
	 	Git commit:	7390fc6
	 	Built:	Tue Feb 27 22:17:40 2018
	 	OS/Arch:	linux/amd64

	Server:
	  Engine:
	  	Version:	17.12.1-ce
	  	API version:	1.35 (minimum version 1.12)
	  	Go version:	go1.9.4
	  	Git commit:	7390fc6
	  	Built:	Tue Feb 27 22:16:13 2018
	  	OS/Arch:	linux/amd64
	  	Experimental:	false
	```


Once Docker is installed, the dockerfiles for each piece of software can be found in the ```./dockerfiles``` folder, and are built using versioning to maintain reproducibility. Docker Images can be obtained locally in 1 of 2 ways:

__Pull__

	```
	make pull
	```

__Build__

	```
	make build
	```	

An example build and deploy of all the images can be found in ./dockerfiles/make.log


## Usage

command Example:

```
nextflow main.nf --sample "single" --reference "hg38" --strand "unstranded" --ercc --fullCov -w "/dump/work" -profiles standard
```

This command will read files from ```./input```, run as __single__ samples, with __hg38__ as the reference, __unstranded__ samples, __ercc__ process will run and so will the __fullCoverage__ process. The work directory has been set to the path ```/dump/work``` and the __standard__ profile has been selected.

This RNA-Seq pipeline can handle several combinations of single/paired read, human/mouse/rat reference, and unstranded'forward/reverse options. The following __mandatory__ options are available:

__sample__

	```
	--sample "single"
	  OR
	--sample "paired"
	```

__reference__

	```
	Human:
	--reference "hg38"
	  OR
	--reference "hg19"

	Mouse:
	--reference "mm10"

	Rat:
	--reference "rn6"
	```

__stranded__

	```
	Stranded:
	--strand "forward"
	  OR
	--strand "reverse"

	Unstranded:
	--strand "unstranded"
	```


### File Naming

NOTICE: File names can not contain "." in the name because the pipeline operates on file names by splitting along the "." to determine prefixes. Change all "." to "\_"

The pipeline can handle merging required samples and files should be named in the following format:

### Merging Required

__single__

	```
	Read 1: "{prefix}_read1.fastq.gz"
	Read 2: "{prefix}_read2.fastq.gz"
	```

__paired__

	```
	Pair 1, Read 1: "{prefix}_1_read1.fastq.gz"
	Pair 1, Read 2: "{prefix}_1_read2.fastq.gz"
	Pair 2, Read 1: "{prefix}_2_read1.fastq.gz"
	Pair 2, Read 2: "{prefix}_2_read2.fastq.gz"
	```


To see what other commands and options the pipeline can handle, type:

```
nextflow main.nf --help
```



## Software

The following software versions are used in this pipeline. For the Ubuntu Base and R Base images, look to the dockerfiles to see specifics on installed packages.

__R Base__

	```
	R_IMAGE = r3.4.3_base
	```

__Ubuntu Base__

	```
	UBUNTU_BASE_IMAGE = ubuntu16.04_base
	```

__Kallisto__

	```
	ERCC_IMAGE = kallisto_v0.43.1
	```

__FastQC__

	```
	QUALITY_IMAGE = fastqc_v0.11.5
	```

__Trimmomatic__

	```
	TRIM_IMAGE = trimmomatic-0.36
	```

__HISAT__

	```
	HISAT_IMAGE = hisat2-2.0.4
	```

__RSeQC__

	```
	RSEQC_IMAGE = rseqc-v2.6.4
	```

__Samtools__

	```
	SAMTOOLS_IMAGE = samtools-1.3.1
	```

__Salmon__

	```
	SALMON_IMAGE = salmon-0.9.1
	```

__Regtools__

	```
	REGTOOLS_IMAGE = regtools-0.3.0
	```

__SubRead__

	```
	SUBREAD_IMAGE = subread-1.6.0
	```

__Wiggletools__

	```
	WIGGLETOOLS_IMAGE = wiggletools-1.2
	```

__BCFTools__

	```
	BCFTOOLS_IMAGE = bcftools-1.3.1
	```

__VCFTOols__

	```
	VCFTOOLS_IMAGE = vcftools-v0.1.15
	```


#### Potential Issues

1) R-Base Dockerfile: This dockerfile holds all of the R packages needed throughout the pipeline. In order to maintain a specific version combination for required software, all packages are manually installed from source. However, this causes the dockerfile to have numerous commit layers, and risks reaching maximum depth ~125 layers. In this case, packages can be consolidated and installed by name (rather than source link) to decrease the number of layers. OR, the existing docker image can be squashed using --squash (with docker in --experimental mode) to flatten the existing image
2) The dockerfiles declare specific versions of the required software, and as software becomes outdated and links need to be updated, the dockerfiles will become less likely to build successfully. In this case its much easier to simply run ```make pull``` in the dockerfiles directory to pull the current working version of each image.

#### Further Development

When the pipeline is completed, the output log files in ```./pipeline_log``` will show the process details for each task. From this information, resource requirements can be defined more accurately to speed up the pipeline execution.