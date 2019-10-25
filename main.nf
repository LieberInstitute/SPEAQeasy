#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
===================================================================================================================================
					LIEBER INSTITUTE JAFFE-LAB	 R N A - S E Q	 A N A L Y S I S	 P I P E L I N E  
===================================================================================================================================
 RNA-Seq Multi-Input Analysis Pipeline. Nextflow Version: Started December 2017.
 #### Homepage / Documentation
 https://github.com/LieberInstitute/RNAsp
 #### Authors
 ##### Original Pipeline
 Emily Burke <emily.burke@libd.org>
 Leonardo Collado-Tores <lcolladotor@gmail.com>
 Andrew Jaffe <andrew.jaffe@libd.org>
 BaDoi Phan <address@email.com>
 ##### Nextflow Version
 Jacob Leonard <leonard.jacob09@gmail.com>
 Israel Aguilar <iaguilaror@gmail.com>
 Violeta Larios <siedracko@gmail.com>
 Nick Eagles <nick.eagles@libd.org>
 ###### References
 https://github.com/SciLifeLab/NGI-smRNAseq
-----------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------
 Pipeline Overview:

   Preprocessing:
	-   Ia: Download Assembly FA
	-   Ib: Build HISAT Index
	-  IIa: Download GENCODE GTF File
	-  IIb: Build Bed File
	- IIIa: Download Salmon TX FA
	- IIIb: Build Salmon Index
	Sample Processing:
	-   A: File Merging (Optional)
	-   B: ERCC Quality Analysis (Optional)
	-   C: Sample Manifest
	-   1: FastQC Quality Analysis
	-  2a: Adaptive Trimming Filter (Sample Dependent)
	-  2b: File Trimming (Sample Dependent)
	-  2c: FastQC Trimmed Quality Analysis (Sample Dependent)
	-  3a: Hisat2 Index Creation
	-  3b: Convert Sam to Bam
	-  3c: Infer Experiment
	-  3d: Infer Strandness
	-  4a: Feature Counts
	-  4b: Primary Alignments
	-  4c: Junctions
	-  5a: Coverage
	-  5b: WigtoBigWig
	-  5c: Mean Coverage
	-   6: Salmon TXQuant
	-  7a: Create Counts Objects
	-  7b: Create Coverage Objects
	-  8a: Call Variants
	-  8b: Merge Called Variants
	-   9: Expressed Regions
-----------------------------------------------------------------------------------------------------------------------------------
*/

def helpMessage() {
	log.info"""
	=============================================================
	 LIEBER INSTITUTE JAFFE-LAB RNA-seq : RNA-Seq Multi-Input Analysis v${version}
	=============================================================
	Usage:
	The typical command for running the pipeline is as follows:

	nextflow run main.nf --sample "single" --strand "unstranded" --reference "hg19" --ercc --fullCov -profile standard

	NOTES: The pipeline accepts single or paired end reads. These reads can be stranded or unstranded.

	NOTE: File names can not have "." in the name due to the prefix functions in between process

  nextflow run main.nf {CORE} {OPTIONS}
    {CORE}:
    --sample "single/paired"
      single <- reads files individually
      paired <- reads files paired
    --reference
      hg38 <- uses human reference hg38
      hg19 <- uses human reference hg19
      mm10 <- uses mouse reference mm10
      rn6  <- uses rat reference rn6
    --strand "forward"/"reverse"/"unstranded"
      forward <- uses forward stranding
      reverse <- uses reverse stranding
      unstranded <- uses pipeine inferencing
    {OPTIONS}:
    --ercc  <- performs ercc quantification
    --fullCov <- performs fullCov R analysis
    --force_trim <- performs trimming on all inputs
  
    --help <- shows this message

	Mandatory Options:
	-----------------------------------------------------------------------------------------------------------------------------------
	--sample		Runs the pipeline for "single" or "paired" end reads
	-----------------------------------------------------------------------------------------------------------------------------------
	--strand		Runs pipeline for "unstranded, forward, reverse" read types
	-----------------------------------------------------------------------------------------------------------------------------------
	--reference		Select the desired reference for the run (hg38, hg19, mm10, rn6)
	-----------------------------------------------------------------------------------------------------------------------------------

	Optional Parameters:
	-----------------------------------------------------------------------------------------------------------------------------------
	--experiment	Name of the experiment being run (ex: "alzheimer"). Defaults to "Jlab_experiment"
	-----------------------------------------------------------------------------------------------------------------------------------
	--prefix		Defines the prefix of the input files (not used to detect files)
	-----------------------------------------------------------------------------------------------------------------------------------
	--input			Defines the input folder for the files. Defaults to "./input" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--output		Defines the output folder for the files. Defaults to "./results" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--unalign		Give the option to not algin the reads against a reference in HISAT step. Defaults to false 
	-----------------------------------------------------------------------------------------------------------------------------------
	--annotation	Path to the folder containing pipeline annotations. Defaults to "./Annotations" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--genotype		Path to the folder containing pipeline genotypes. Defaults to "./Genotyping" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--ercc			Flag to enable ERCC quantification with Kallisto
	-----------------------------------------------------------------------------------------------------------------------------------
	--fullCov		Flag to perform full coverage in step 7b
	-----------------------------------------------------------------------------------------------------------------------------------
	--small_test	Runs the pipeline as a small test run on sample files located in the test folder
	-----------------------------------------------------------------------------------------------------------------------------------
    --force_trim    Include to perform trimming on all inputs, rather than just those failing QC due to detected adapter content
    -----------------------------------------------------------------------------------------------------------------------------------
	""".stripIndent()
}
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = "0.8.0"

// Show help message
params.help = false
if (params.help){
	helpMessage()
	exit 0
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define default values for params
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

params.experiment = "Jlab_experiment"
params.prefix = "prefix"
params.sample = ""
params.strand = ""
params.unalign = false
params.reference = ""
params.annotation = "${workflow.projectDir}/Annotation"
params.genotype = "${workflow.projectDir}/Genotyping"
params.output = "${workflow.projectDir}/results"
params.scripts = "${workflow.projectDir}/scripts"
params.ercc = false
params.fullCov = false
params.small_test = false
params.force_trim = false
workflow.runName = "RNAsp_run"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Validate Inputs
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Sample Selection Validation
if (params.sample != "single" && params.sample != "paired") {
	exit 1, "Sample Type Not Provided or Invalid Choice. Please Provide a Valid Sample Type. Valid types are single or paired "
}

// Strand Selection Validation
if (params.strand != "forward" && params.strand != "reverse" && params.strand != "unstranded") {
	exit 1, "Strand Type Not Provided or Invalid Choice. Please Provide a Valid Strand Type. Valid types are unstranded, forward or reverse"
}

// Reference Selection Validation
if (params.reference != "hg19" && params.reference != "hg38" && params.reference != "mm10" && params.reference != "rn6") {
	exit 1, "Error: enter hg19 or hg38, mm10 for mouse, or rn6 for rat as the reference."
}

// Get species name from genome build name
if (params.reference == "hg19" || params.reference == "hg38") {
    params.reference_type = "human"
} else if (params.reference == "mm10") {
    params.reference_type = "mouse"
} else {
    params.reference_type = "rat"
}

//  Path to small test files
if (params.small_test) {
  if (params.strand == "unstranded") {
    params.input = "${workflow.projectDir}/test/$params.reference_type/${params.sample}/unstranded"
  } else {
    params.input = "${workflow.projectDir}/test/$params.reference_type/${params.sample}/stranded"
  }
} else {
  params.input = "${workflow.projectDir}/input"
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Strand Option Parameters
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  Single vs. paired-end-related command-line options
if (params.sample == "single") {
    params.trim_sample = "SE"
    params.kallisto_strand = "--single -l ${params.kallisto_len_mean} -s ${params.kallisto_len_sd}"
} else {
    params.trim_sample = "PE"
    params.kallisto_strand = ""
}

//  Strandness-related command-line options
if (params.strand == "unstranded") {
    params.feature_strand = "0"
    params.hisat_strand = ""
} else if (params.strand == "forward") {
    params.feature_strand = "1"
    params.kallisto_strand += " --fr-stranded"
    if (params.sample == "single") {
        params.hisat_strand = "--rna-strandness F"
    } else {
        params.hisat_strand = "--rna-strandness FR"
    }
} else {
    params.feature_strand = "2"
    params.kallisto_strand += " --rf-stranded"
    if (params.sample == "single") {
        params.hisat_strand = "--rna-strandness R"
    } else {
        params.hisat_strand = "--rna-strandness RF"
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define External Scripts
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

infer_strandness_script = file("${params.scripts}/step3b_infer_strandness.R")
prep_bed_script = file("${params.scripts}/prep_bed.R")
bed_to_juncs_script = file("${params.scripts}/bed_to_juncs.py")
expressed_regions_script = file("${params.scripts}/step9-find_expressed_regions.R")
fullCov_script = file("${params.scripts}/create_fullCov_object.R")

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define Reference Paths/Scripts + Reference Dependent Parameters
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ERCC
if (params.ercc) {
  erccidx = file("${params.annotation}/ERCC/ERCC92.idx")
  ercc_actual_conc = file("${params.annotation}/ercc_actual_conc.txt")
}

if (params.reference == "hg38") {
	
	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz"
	params.fa = "GRCh38.primary_assembly.genome.fa"
	params.hisat_prefix = "hisat2_GRCh38primary"
	params.hisat_assembly = "GENCODE/GRCh38_hg38/assembly"

	// Step 4: gencode gtf
	params.gencode_gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"
	params.gencode_gtf = "gencode.v25.annotation.gtf"
	params.feature_output_prefix = "Gencode.v25.hg38"

	// Step 5: python coverage
    //##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
	chr_sizes = file("${params.annotation}/chrom_sizes/hg38.chrom.sizes.gencode")

	// Step 6: salmon
    //##TODO(iaguilar): Explain why step 6 is enabled if reference is hg38...  (Doc ######)
	params.step6 = true
	params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.transcripts.fa.gz"
	params.tx_fa = "gencode.v25.transcripts.fa"
	params.salmon_prefix = "salmon_0.8.2_index_gencode.v25.transcripts"
	params.salmon_assembly = "GENCODE/GRCh38_hg38/transcripts"

	// Step 7: Make R objects	
	junction_annotation_gencode = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_hg38_gencode_v25.rda")
	junction_annotation_ensembl = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_hg38_ensembl_v85.rda")
	junction_annotation_genes = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_hg38_refseq_grch38.rda")
	feature_to_tx_gencode = Channel.fromPath("${params.annotation}/junction_txdb/feature_to_Tx_hg38_gencode_v25.rda")
	feature_to_tx_ensembl = Channel.fromPath("${params.annotation}/junction_txdb/feature_to_Tx_ensembl_v85.rda")
    exon_maps_by_coord_hg38 = Channel.fromPath("${params.annotation}/junction_txdb/exonMaps_by_coord_hg38_gencode_v25.rda")
	create_counts = file("${params.scripts}/create_count_objects-human.R")

	// Step 8: call variants
    //##TODO(iaguilar): Explain why step 8 is enabled if reference is hg38...  (Doc ######)
	params.step8 = true
	snvbed = Channel.fromPath("${params.genotype}/common_missense_SNVs_hg38.bed")

}
if (params.reference == "hg19") {
	
	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
	params.fa = "GRCh37.primary_assembly.genome.fa"
	params.hisat_prefix = "hisat2_GRCh37primary"
	params.hisat_assembly = "GENCODE/GRCh37_hg19/assembly"
	
	// Step 4: gencode gtf
	params.gencode_gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz"
	params.gencode_gtf = "gencode.v25lift37.annotation.gtf"
	params.feature_output_prefix = "Gencode.v25lift37.hg19"

	// Step 5: python coverage
    //##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
	chr_sizes = file("${params.annotation}/chrom_sizes/hg19.chrom.sizes.gencode")

	// Step 6: salmon
    //##TODO(iaguilar): Explain why step 6 is enabled if reference is hg19...  (Doc ######)
	params.step6 = true
	params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.transcripts.fa.gz"
	params.tx_fa = "gencode.v25lift37.transcripts.fa"
	params.salmon_prefix = "salmon_0.8.2_index_gencode.v25lift37.transcripts"
	params.salmon_assembly = "GENCODE/GRCh37_hg19/transcripts"

	// Step 7: Make R objects
	junction_annotation_gencode = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_hg19_gencode_v25lift37.rda")
	junction_annotation_ensembl = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_hg19_ensembl_v75.rda")
	junction_annotation_genes = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_hg19_refseq_grch37.rda")
	feature_to_tx_gencode = Channel.fromPath("${params.annotation}/junction_txdb/feature_to_Tx_hg19_gencode_v25lift37.rda")
	feature_to_tx_ensembl = Channel.fromPath("${params.annotation}/junction_txdb/feature_to_Tx_ensembl_v75.rda")
	create_counts = file("${params.scripts}/create_count_objects-human.R")

	// Step 8: call variants
    //##TODO(iaguilar): Explain why step 8 is enabled if reference is hg19...  (Doc ######)
	params.step8 = true
	snvbed = Channel.fromPath("${params.genotype}/common_missense_SNVs_hg19.bed")

}
if (params.reference == "mm10") {

	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/GRCm38.primary_assembly.genome.fa.gz"
	params.fa = "GRCm38.primary_assembly.genome.fa"
	params.hisat_prefix = "GRCm38_mmhisat2_GRCm38primary"
	params.hisat_assembly = "GENCODE/GRCm38_mm10/assembly"

	// Step 4: gencode gtf
	params.gencode_gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/gencode.vM11.annotation.gtf.gz"
	params.gencode_gtf = "gencode.vM11.annotation.gtf"
	params.feature_output_prefix = "Gencode.M11.mm10"

	// Step 5: python coverage
    //##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
	chr_sizes = file("${params.annotation}/chrom_sizes/mm10.chrom.sizes.gencode")

	// Step 6: salmon
    //##TODO(iaguilar): Explain why step 6 is enabled if reference is mm10...  (Doc ######)
	params.step6 = true
	params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/gencode.vM11.transcripts.fa.gz"
	params.tx_fa = "gencode.vM11.transcripts.fa"
	params.salmon_prefix = "salmon_0.8.2_index_gencode.vM11.transcripts"
	params.salmon_assembly = "GENCODE/GRCm38_mm10/transcripts"

	// Step 7: Make R objects
	junction_annotation_gencode = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_mm10_gencode_vM11.rda")
	junction_annotation_ensembl = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_mm10_ensembl_v86.rda")
	create_counts = file("${params.scripts}/create_count_objects-mouse.R")

    // Step 8: call variants
    //##TODO(iaguilar): Explain why step 8 is disabled if reference is mm10...  (Doc ######)
	params.step8 = false
}
if (params.reference == "rn6") {

	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ensembl.org/pub/release-86/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
	params.fa = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
	params.hisat_prefix = "hisat2_Rnor6.0toplevel"
	params.hisat_assembly = "ensembl/Rnor_6.0"

	// Step 4: gencode gtf (ensembl for rn6)
	params.gencode_gtf_link = "ftp://ftp.ensembl.org/pub/release-86/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.86.gtf.gz"
	params.gencode_gtf = "Rattus_norvegicus.Rnor_6.0.86.gtf"
	params.feature_output_prefix = "Rnor_6.0.86"

	// Step 5: python coverage
    //##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
	chr_sizes = file("${params.annotation}/chrom_sizes/rn6.chrom.sizes.ensembl")
	
	// Step 6: Salmon
    //##TODO(iaguilar): Explain why step 6 is enabled if reference is rn6...  (Doc ######)
	params.step6 = false

	// Step 7: Make R objects
	junction_annotation_ensembl = Channel.fromPath("${params.annotation}/junction_txdb/junction_annotation_rn6_ensembl_v86.rda")
	create_counts = file("${params.scripts}/create_count_objects-rat.R")

	//Step 8: call variants
    //##TODO(iaguilar): Explain why step 8 is enabled if reference is rn6...  (Doc ######)
	params.step8 = false

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Utilities for retrieving info from filenames
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_prefix(f) {
  //  Remove these regardless of position in the string
  String blackListAny = "_summary|_TR|_TNR|_trimmed|_reverse_paired|_1|_2|_reverse_unpaired|_forward_paired|_forward_unpaired|_hisat_out"
  
  f.name.toString()
   .replaceAll("_read1.", ".")
   .replaceAll("_read2.", ".")
   .replaceAll(".Forward", "_Forward")
   .replaceAll(".Reverse", "_Reverse")
   .tokenize('.')[0]
   .replaceAll(blackListAny, "")
}

def get_read_type(f) {
  baseName = f.name.toString().tokenize('.')[0]
  if (baseName.endsWith('_Forward')) {
    return('Forward')
  } else if (baseName.endsWith('_Reverse')) {
    return('Reverse')
  } else {
    return('Unstranded')
  }
}

def get_file_ext(f) {
  if (f.name.toString().tokenize(".")[-1] == "gz") {
    return('.fastq.gz')
  } else {
    return('.fastq')
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Summary of Defined Variables
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Header log info
log.info "============================================================="
log.info " LIBD-RNAseq : Multi-Input RNA-Seq Best Practice v${version}"
log.info "============================================================="
def summary = [:]
summary['Run Name']			= workflow.runName
summary['Sample']			  = params.sample
summary['Reference']		   = params.reference
summary['Strand']			  = params.strand
summary['Annotation']		 = params.annotation
summary['Genotype']		   = params.genotype
summary['Input']			   = params.input
summary['scripts']      = params.scripts
if(params.ercc) summary['ERCC Index'] = erccidx
if(params.experiment) summary['Experiment'] = params.experiment
if(params.unalign) summary['Align'] = "True"
if(params.fullCov) summary['Full Coverage'] = "True"
summary['Small test selected'] = params.small_test
summary['Output dir']		  = params.output
summary['Working dir']		 = workflow.workDir
summary['Current home']		= "$HOME"
summary['Current user']		= "$USER"
summary['Current path']		= "$PWD"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================="

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN PIPELINE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* ###################################

	STARTS REFERENCES BUILDING BLOCK 

################################### */

/*
 * Step Ia: GENCODE Assembly FA
 */

// Define the output directory path for Genome Reference dir, and hisat2 index dir, based on previously defined parameters
params.hisat_idx_output = "${params.annotation}/${params.hisat_assembly}"

// Uses "storeDir" to download assembly file if required, and otherwise output the cached file
// if it already exists
process pullGENCODEassemblyfa {

  tag "Downloading Assembly FA File: ${params.fa}"
  storeDir "${params.hisat_idx_output}/fa"

	output:
	  file "${params.fa}" into reference_assembly, variant_assembly

	shell:
	  '''
    wget "!{params.fa_link}"
    gunzip "!{params.fa}.gz"
    '''
}

/*
 * Step Ib: Build HISAT Index
 */

// Uses "storeDir" to build HISAT2 index only when it doesn't exist, and output the cached
// files if they do already exist
process buildHISATindex {
		
  tag "Building HISAT2 Index: ${params.hisat_prefix}"
  storeDir "${params.hisat_idx_output}/index"

  input:
    file reference_fasta from reference_assembly

  output:
    file("${params.hisat_prefix}.*") into hisat_index_built

  shell:
    '''
    !{params.hisat2build} -p !{task.cpus} !{reference_fasta} !{params.hisat_prefix}
    '''
}

// Channel post-processing
hisat_index_built // get every *.ht2 file in this channel
	.toList() // group every *.ht2 item into a single list
	.set{ hisat_index } // pass *.ht2 list to a new channel

/*
 * Step IIa: GENCODE GTF Download
 */

// Define the output directory path for GTF file, and the prep bed, based on previously defined parameters
params.gencode_gtf_out = "${params.annotation}/RSeQC/${params.reference}"

// Uses "storeDir" to download gtf only when it doesn't exist, and output the cached
// file if it does already exist
process pullGENCODEgtf {

  tag "Downloading GTF File: ${params.gencode_gtf}"
  storeDir "${params.gencode_gtf_out}/gtf"

  output:
    file "${params.gencode_gtf}" into gencode_gtf, create_counts_gtf, gencode_feature_gtf

  shell:
    '''
    wget "!{params.gencode_gtf_link}"
    gunzip "!{params.gencode_gtf}.gz"
    '''
}

/*
 * Step IIb: Build Bed File
 */

// Uses "storeDir" to build bed file only when it doesn't exist, and output the cached
// file if it does already exist
process buildPrepBED {
	
  tag "Building Bed File: ${params.reference}"
  storeDir "${params.gencode_gtf_out}/bed"

  input:
    file gencode_gtf from gencode_gtf
    file prep_bed_script from prep_bed_script

  output:
    file("${params.reference}.bed") into bedfile

  shell:
    '''
    !{params.Rscript} !{prep_bed_script} -f !{gencode_gtf} -n !{params.reference}
    '''
}

// for hg38, hg19, and mm10, step 6 is enabled by params.step6 = true 
// during Define Reference Paths/Scripts + Reference Dependent Parameters
if (params.step6) {

	params.salmon_idx_output = "${params.annotation}/${params.salmon_assembly}"

	/*
	 * Step IIIa: GENCODE TX FA Download
	 */
		
    // Uses "storeDir" to download files only when they don't exist, and output the cached
    // files if they do already exist
    process pullGENCODEtranscripts {
			
      tag "Downloading TX FA File: ${params.tx_fa}"
      storeDir "${params.salmon_idx_output}/fa"

      output:
        file("${params.tx_fa}") into transcript_fa

      shell:
        '''
        wget !{params.tx_fa_link}
        gunzip !{params.tx_fa}.gz
        '''
    }

	/*
	 * Step IIIb: Salmon Transcript Build
	 */

    // Uses "storeDir" to build salmon index only if the pre-built file is not present; outputs
    // this cached file otherwise
	process buildSALMONindex {

    tag "Building Salmon Index: ${params.salmon_prefix}"
    storeDir "${params.salmon_idx_output}/salmon"

    input:
      file tx_file from transcript_fa

    output:
      file("${params.salmon_prefix}") into salmon_index_built

    script:
      """
      ${params.salmon} index -t $tx_file -i ${params.salmon_prefix} -p $task.cpus --type quasi -k ${params.salmon_min_read_len}
      """
    }

    // Post-processing of built index for use as input to TXQuant process
    salmon_index_built
      .collect()
      .set{ salmon_index }
}

//  Step A: Merge files as necessary

samples_manifest = file("${params.input}/samples.manifest")
merge_script = file("${params.scripts}/step00-merge.R")
process Merging {
  tag "Performing merging if/where necessary"
  
  input:
    file original_manifest from samples_manifest
    file merge_script from merge_script
  output:
    file "out/*" into merged_inputs_flat
  
  shell:
    '''
    !{params.Rscript} !{merge_script} -s !{original_manifest} -o $(pwd)/out -c !{task.cpus}
    '''
}

//  Group both reads together for each sample, if paired-end, and assign each sample a prefix
if (params.sample == "single") {

    merged_inputs_flat
        .flatten()
        .map{file -> tuple(get_prefix(file), file) }
        .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
        .set{ temp_inputs }
} else {

    merged_inputs_flat
        .toSortedList()
        .flatten()
        .map{file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
        .set{ temp_inputs }
}

//  Copy the processed input channel to channels used by dependent processes
if (params.ercc) { 
  temp_inputs.into{ fastqc_untrimmed_inputs; adaptive_trimming_fastqs; manifest_creation; salmon_inputs; ercc_inputs }
} else {
  temp_inputs.into{ fastqc_untrimmed_inputs; adaptive_trimming_fastqs; manifest_creation; salmon_inputs }
}

/*
 * Step B: Run the ERCC process if the --ercc flag is specified
 */
 
if (params.ercc) {

  process ERCC {
		
    tag "Prefix: $ercc_prefix"
    publishDir "${params.output}/ercc/${ercc_prefix}",'mode':'copy'

    input:
      file erccidx from erccidx
      set val(ercc_prefix), file(ercc_input) from ercc_inputs

    output:
      file("${ercc_prefix}_abundance.tsv") into ercc_abundances

    script:
      """
      ${params.kallisto} quant -i $erccidx -t $task.cpus -o . ${params.kallisto_strand} $ercc_input \
      && cp abundance.tsv ${ercc_prefix}_abundance.tsv
      """
  }
}

/*
 * Step C: Sample Manifest
 */

man_info_script = file("${params.scripts}/find_sample_info.R")
process Manifest {
	
	tag "Validating manifest and accounting for merged files"
	publishDir "${params.output}/manifest", mode:'copy'

	input:
	  file original_manifest from samples_manifest
      file man_info_script from man_info_script

	output:
	  file "samples.manifest" into counts_samples_manifest, fullCov_samples_manifest

	shell:
	  '''
      !{params.Rscript} !{man_info_script} -s !{original_manifest} -o "."
	  '''
}

/*
 * Step 1: Untrimmed Quality Report
 */

process QualityUntrimmed {

	tag "Prefix: $untrimmed_prefix"
	publishDir "${params.output}/FastQC/Untrimmed",mode:'copy'

	input:
	set val(untrimmed_prefix), file(fastqc_untrimmed_input) from fastqc_untrimmed_inputs 

	output:
	file "*"
	file "*_summary.txt" into quality_reports, count_objects_quality_reports
	file "*_fastqc_data.txt" into count_objects_quality_metrics

	script:
	if (params.sample == "single") {
		copy_command = "cp ${untrimmed_prefix}_fastqc/summary.txt ${untrimmed_prefix}_summary.txt"
		data_command = "cp ${untrimmed_prefix}_fastqc/fastqc_data.txt ${untrimmed_prefix}_fastqc_data.txt"
	} else {
		copy_command = "cp ${untrimmed_prefix}_1_fastqc/summary.txt ${untrimmed_prefix}_1_summary.txt && cp ${untrimmed_prefix}_2_fastqc/summary.txt ${untrimmed_prefix}_2_summary.txt"
		data_command = "cp ${untrimmed_prefix}_1_fastqc/fastqc_data.txt ${untrimmed_prefix}_1_fastqc_data.txt && cp ${untrimmed_prefix}_2_fastqc/fastqc_data.txt ${untrimmed_prefix}_2_fastqc_data.txt"
	}
	"""
	${params.fastqc} $fastqc_untrimmed_input --extract
	$copy_command
	$data_command
	"""
}

/*
 * Step 2a: Adaptive Trimming
 */

if (params.sample == "single") {

	quality_reports
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .join(adaptive_trimming_fastqs)
	  .ifEmpty{ error "Cannot Find Combined Quality and Trimming Channel for Single Adaptive Trimming" }
	  .set{ adaptive_trimming_single_inputs }

	process AdaptiveTrimSingleReads {
	  
	  tag "Prefix: $single_adaptive_prefix"
	  publishDir "${params.output}/Adaptive_Trim",'mode':'copy'

	  input:
	  set val(single_adaptive_prefix), file(single_adaptive_summary), file(single_adaptive_fastq) from adaptive_trimming_single_inputs

	  output:
	  file "*" into trimming_fastqs, no_trimming_fastqs

	  shell:
	  single_quality_report = single_adaptive_prefix.toString() + "_summary.txt"
	  single_trimming_input = single_adaptive_prefix.toString() + ".f*q*"
      suffix = get_file_ext(single_adaptive_fastq)
	  '''
	  export result=$(grep "Adapter Content" !{single_quality_report} | cut -f1)
	  if [ $result == "FAIL" ] || [ !{params.force_trim} == "true" ] ; then
		  mv !{single_trimming_input} "!{single_adaptive_prefix}_TR!{suffix}"
	  else
		  mv !{single_trimming_input} "!{single_adaptive_prefix}_TNR!{suffix}"
	  fi
	  '''
	}

} else { // paired

	quality_reports
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .groupTuple()
	  .join(adaptive_trimming_fastqs)
	  .ifEmpty{ error "Cannot Find Combined Quality and Trimming Channel for Paired Adaptive Trimming" }
	  .set{ adaptive_trimming_paired_inputs }

	process AdaptiveTrimPairedReads {
	  
	  tag "Prefix: $paired_adaptive_prefix"
	  publishDir "${params.output}/Adaptive_Trim",mode:'copy'

	  input:
	  set val(paired_adaptive_prefix), file(paired_adaptive_summary), file(paired_adaptive_fastq) from adaptive_trimming_paired_inputs

	  output:
	  file "*" into trimming_fastqs, no_trimming_fastqs

	  shell:
      suffix = get_file_ext(paired_adaptive_fastq[0])
      prefix_str = paired_adaptive_prefix.toString()
	  '''
	  export result1=$(grep "Adapter Content" !{prefix_str}_1_summary.txt | cut -c1-4)
	  export result2=$(grep "Adapter Content" !{prefix_str}_2_summary.txt | cut -c1-4)
	  if [ $result1 == "FAIL" ] || [ $result2 == "FAIL" ] || [ !{params.force_trim} == "true" ]; then
		  cp !{prefix_str}_1.f*q* "!{prefix_str}_1_TR!{suffix}"
		  cp !{prefix_str}_2.f*q* "!{prefix_str}_2_TR!{suffix}"
	  else
		  cp !{prefix_str}_1.f*q* "!{prefix_str}_1_TNR!{suffix}"
		  cp !{prefix_str}_2.f*q* "!{prefix_str}_2_TNR!{suffix}"
	  fi
	  '''
	}
}

/*
 * Modify the Trimming Input Channel 
 */

if (params.sample == "single") {

	trimming_fastqs
	  .flatten()
	  .filter{ file -> file.name.toString() =~ /_TR.*/ }
	  .map{ file -> tuple(get_prefix(file), file) }
	  .set{ trimming_inputs }

	no_trimming_fastqs
	  .flatten()
	  .filter{ file -> file.name.toString() =~ /_TNR.*/ }
	  .map{ file -> tuple(get_prefix(file), file) }
	  .set{ no_trim_fastqs }

} else { // paired

	trimming_fastqs
	  .flatten()
	  .filter{ file -> file.name.toString() =~ /_TR.*/ }
	  .toSortedList()
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .groupTuple()
	  .set{ trimming_inputs }

	no_trimming_fastqs
	  .flatten()
	  .filter{ file -> file.name.toString() =~ /_TNR.*/ }
	  .toSortedList()
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .groupTuple()
	  .set{ no_trim_fastqs }
}

/*
 * Step 2b: Trimming 
 */

process Trimming {
	
	tag "Prefix: $trimming_prefix"
	publishDir "${params.output}/trimmed_fq",'mode':'copy'

	input:
	set val(trimming_prefix), file(trimming_input) from trimming_inputs

	output:
	file "*.f*q*" into trimmed_fastqc_inputs, trimmed_hisat_inputs

	script:
	sample_option = "${params.trim_sample}"
	if (params.sample == "single") {
		output_option = "${trimming_prefix}_trimmed.f*q*"
	}
	if (params.sample == "paired") {
		output_option = "${trimming_prefix}_trimmed_forward_paired.f*q* ${trimming_prefix}_trimmed_forward_unpaired.f*q* ${trimming_prefix}_trimmed_reverse_paired.f*q* ${trimming_prefix}_trimmed_reverse_unpaired.f*q*"
	}
	"""
    #  This solves the problem of trimmomatic and the adapter fasta
    #  needing hard paths, even when on the PATH.
    if [ ${params.use_long_paths} == "true" ]; then
        trim_jar=${params.trimmomatic}
        adapter_fa=${params.adapter_fasta}
    else
        trim_jar=\$(which ${params.trimmomatic})
        adapter_fa=\$(which ${params.adapter_fasta})
    fi

	java -Xmx512M \
	-jar \$trim_jar \
	$sample_option \
	-threads $task.cpus \
	-phred33 \
	$trimming_input \
	$output_option \
	ILLUMINACLIP:\$adapter_fa:2:30:10:1 \
	LEADING:3 \
	TRAILING:3 \
	SLIDINGWINDOW:4:15 \
	MINLEN:75
	"""
}


/*
 * Step 2c: Run FastQC Quality Check on Trimmed Files
 */

process QualityTrimmed {

	tag "Prefix: $fastq_name"
	publishDir "${params.output}/FastQC/Trimmed",'mode':'copy'

	input:
	file fastqc_trimmed_input from trimmed_fastqc_inputs

	output:
	file "*"

	script:
    fastq_name = get_prefix(fastqc_trimmed_input)
	"""
	$params.fastqc $fastqc_trimmed_input --extract
	"""
}

/*
 * Step 3a: Hisat Sam File
 */

if (params.sample == "single") {

    //Here trimmed and not trimmed data is mixed in a channel to ensure the flow of the pipeline
	trimmed_hisat_inputs
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .mix(no_trim_fastqs)
	  .ifEmpty{ error "Single End Channel for HISAT is empty" }
	  .set{ single_hisat_inputs }

  process SingleEndHISAT {

	  tag "Prefix: $single_hisat_prefix"
	  publishDir "${params.output}/HISAT2_out",mode:'copy'

	  input:
	  file hisat_index from hisat_index
	  set val(single_hisat_prefix), file(single_hisat_input) from single_hisat_inputs

	  output:
	  file "*_hisat_out.sam" into hisat_single_output
	  file "*"
	  file "*_align_summary.txt" into alignment_summaries

	  shell:
	  hisat_full_prefix = "${params.annotation}/${params.hisat_assembly}/index/${params.hisat_prefix}"
	  strand = "${params.hisat_strand}"
	  // Phred Quality is hardcoded, if it will be so, it should be pointed in the README.md
	  '''
	  !{params.hisat2} -p !{task.cpus} -x !{hisat_full_prefix} -U !{single_hisat_input} -S !{single_hisat_prefix}_hisat_out.sam !{strand} --phred33 2> !{single_hisat_prefix}_align_summary.txt
	  '''
	}

	hisat_single_output
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .set{ sam_to_bam_inputs }
     
} else { // sample is paired-end

	no_trim_fastqs
	  .set{ notrim_paired_hisat_inputs }

	process PairedEndNoTrimHISAT {
  
	  tag "Prefix: $paired_notrim_hisat_prefix"
	  publishDir "${params.output}/HISAT2_out",'mode':'copy'

	  input:
	  file hisatidx from hisat_index
	  set val(paired_notrim_hisat_prefix), file(paired_no_trim_hisat) from notrim_paired_hisat_inputs

	  output:
	  file "*_hisat_out.sam" into hisat_paired_notrim_output
	  file "*"
	  file "*_align_summary.txt" into paired_notrim_alignment_summaries

	  shell:
	  hisat_full_prefix = "${params.annotation}/${params.hisat_assembly}/index/${params.hisat_prefix}"
	  strand = "${params.hisat_strand}"
	  sample_1_hisat = paired_notrim_hisat_prefix.toString() + "_1_TNR.f*q*"
	  sample_2_hisat = paired_notrim_hisat_prefix.toString() + "_2_TNR.f*q*"
	  if (params.unalign) {
		  unaligned_opt = "--un-conc ${paired_notrim_hisat_prefix}.fastq"
	  }
	  if (!params.unalign) {
		  unaligned_opt = ""
	  }
	  '''
	  !{params.hisat2} \
	  -p !{task.cpus} \
	  -x !{hisat_full_prefix} \
	  -1 !{sample_1_hisat} \
	  -2 !{sample_2_hisat} \
	  -S !{paired_notrim_hisat_prefix}_hisat_out.sam !{strand} --phred33 \
	  !{unaligned_opt} \
	  2> !{paired_notrim_hisat_prefix}_align_summary.txt
	  '''
  }

	 trimmed_hisat_inputs
	  .flatten()
	  .map{ file -> tuple(get_prefix(file), file) }
	  .groupTuple()
	  .set{ trim_paired_hisat_inputs }

  process PairedEndTrimmedHISAT {
 
	  tag "Prefix: $paired_trimmed_prefix"
	  publishDir "${params.output}/HISAT2_out",'mode':'copy'

	  input:
	  file hisatidx from hisat_index
	  set val(paired_trimmed_prefix), file(paired_trimmed_fastqs) from trim_paired_hisat_inputs

	  output:
	  file "*_hisat_out.sam" into hisat_paired_trim_output
	  file "*"
	  file "*_align_summary.txt" into paired_trim_alignment_summaries

	  shell:
	  hisat_full_prefix = "${params.annotation}/${params.hisat_assembly}/index/${params.hisat_prefix}"
	  strand = "${params.hisat_strand}"
	  forward_paired = paired_trimmed_prefix.toString() + "_trimmed_forward_paired.f*q*"
	  reverse_paired = paired_trimmed_prefix.toString() + "_trimmed_reverse_paired.f*q*"
	  forward_unpaired = paired_trimmed_prefix.toString() + "_trimmed_forward_unpaired.f*q*"
	  reverse_unpaired = paired_trimmed_prefix.toString() + "_trimmed_reverse_unpaired.f*q*"
	  if (params.unalign) {
		  unaligned_opt = "--un-conc ${paired_trimmed_prefix}.fastq"
	  }
	  if (!params.unalign) {
		  unaligned_opt = ""
	  }
	  '''
	  !{params.hisat2} \
	  -p !{task.cpus} \
	  -x !{hisat_full_prefix} \
	  -1 !{forward_paired} \
	  -2 !{reverse_paired} \
	  -U !{forward_unpaired},!{reverse_unpaired} \
	  -S !{paired_trimmed_prefix}_hisat_out.sam !{strand} --phred33 \
	  !{unaligned_opt} \
	  2> !{paired_trimmed_prefix}_align_summary.txt
	  '''
  }
  
  hisat_paired_notrim_output
    .mix(hisat_paired_trim_output)
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .set{ sam_to_bam_inputs }

	paired_trim_alignment_summaries
	  .mix(paired_notrim_alignment_summaries)
	  .flatten()
	  .set{ alignment_summaries } // think this is causing conflicts...
}

/*
 * Step 3b: Sam to Bam 
 */

process SamtoBam {
	
	tag "Prefix: $sam_to_bam_prefix"
	publishDir "${params.output}/HISAT2_out/sam_to_bam",'mode':'copy'

	input:
	set val(sam_to_bam_prefix), file(sam_to_bam_input) from sam_to_bam_inputs

	output:
	set val("${sam_to_bam_prefix}"), file("${sam_to_bam_prefix}*.sorted.bam"), file("${sam_to_bam_prefix}*.sorted.bam.bai") into infer_experiment_inputs, feature_bam_inputs, alignment_bam_inputs, coverage_bam_inputs, full_coverage_bams, count_objects_bam_files, variant_calls_bam

	script:
	original_bam = "${sam_to_bam_prefix}_accepted_hits.bam"
	sorted_bam = "${sam_to_bam_prefix}_accepted_hits.sorted"
	"""
	${params.samtools} view -bh -F 4 $sam_to_bam_input > $original_bam
	${params.samtools} sort -T temporary -@ $task.cpus $original_bam -o ${sorted_bam}.bam
	${params.samtools} index ${sorted_bam}.bam
	"""
}

infer_experiment_inputs
  .combine(bedfile)
  .set{ infer_experiments }

/*
 * Step 3c: Infer Experiment
 */

process InferExperiment {

	tag "Prefix: $infer_prefix"
	publishDir "${params.output}/HISAT2_out/infer_experiment",'mode':'copy'

	input:
	set val(infer_prefix), file(bam_file), file(bam_index), file(bed_file) from infer_experiments

	output:
	file "*"
	file "*_experiment.txt" into infer_experiment_outputs

	shell:
	'''
	!{params.python} $(which !{params.infer_experiment}) \
	-i !{bam_file} \
	-r !{bed_file} \
	1> !{infer_prefix}_experiment.txt \
	2> !{infer_prefix}_experiment_summary_out.txt
	'''
}

infer_experiment_outputs
  .flatten()
  .toSortedList()
  .set{ infer_experiment_output }

/*
 * Step 3d: Infer Strandness
 */

process InferStrandness {

	publishDir "${params.output}/HISAT2_out/infer_strandness/",'mode':'copy'

	input:
	file infer_strandness_script from infer_strandness_script
	file infer_experiment_files from infer_experiment_output

	output:
	file "*"
	file "inferred_strandness_pattern.txt" into inferred_strand_coverage, inferred_strand_mean_coverage, inferred_strand_objects

	shell:
	'''
	!{params.Rscript} !{infer_strandness_script} -p inferred_strandness_pattern.txt
	'''
}

feature_bam_inputs
  .combine(gencode_feature_gtf)
  .set{ feature_counts_inputs }

/*
 * Step 4a: Feature Counts
 */

process FeatureCounts {

	tag "Prefix: $feature_prefix"
	publishDir "${params.output}/Counts",'mode':'copy'

	input:
	set val(feature_prefix), file(feature_bam), file(feature_index), file(gencode_gtf_feature) from feature_counts_inputs

	output:
	file "*"
	file "*.counts*" into sample_counts

	script:
	if (params.sample == "single") {
		sample_option = ""
	} else {
		sample_option = "-p"
	}
	feature_out = "${feature_prefix}_${params.feature_output_prefix}"
	"""
	${params.featureCounts} \
	-s ${params.feature_strand} \
	$sample_option \
	-T $task.cpus \
	-a $gencode_gtf_feature \
	-o ${feature_out}_Genes.counts \
	$feature_bam

	${params.featureCounts} \
	-s ${params.feature_strand} \
	$sample_option \
	-O \
	-f \
	-T $task.cpus \
	-a $gencode_gtf_feature \
	-o ${feature_out}_Exons.counts \
	$feature_bam
	"""
}

/*
 * Step 4b: Primary Alignments
 */

process PrimaryAlignments {

	tag "Prefix: $alignment_prefix"
	publishDir "${params.output}/Counts/junction/primary_aligments",'mode':'copy'

	input:
	set val(alignment_prefix), file(alignment_bam), file(alignment_index) from alignment_bam_inputs

	output:
	set val("${alignment_prefix}"), file("${alignment_prefix}.bam"), file("${alignment_prefix}.bam.bai") into primary_alignments

	script:
	"""
	${params.samtools} view -@ $task.cpus -bh -F 0x100 $alignment_bam > ${alignment_prefix}.bam
	${params.samtools} index ${alignment_prefix}.bam
	"""
}

/*
 * Step 4c: Junctions
 */

process Junctions {

	tag "Prefix: $junction_prefix"
	publishDir "${params.output}/Counts/junction",'mode':'copy'

	input:
	file bed_to_juncs_script from bed_to_juncs_script
	set val(junction_prefix), file(alignment_bam), file(alignment_index) from primary_alignments

	output:
	file "*"
	file("*.count") into regtools_outputs
	//needs to pass the count files to a channel

	shell:
	outjxn = "${junction_prefix}_junctions_primaryOnly_regtools.bed"
	outcount = "${junction_prefix}_junctions_primaryOnly_regtools.count"
	'''
	!{params.regtools} junctions extract -i !{params.juncts_min_intron_len} -o !{outjxn} !{alignment_bam}
	!{params.python} !{bed_to_juncs_script} < !{outjxn} > !{outcount}
	'''
}

/*
 * Step 5a: Coverage
 */

//
if (params.strand == "unstranded") {
	params.strandprefix="*"
} else if (params.strand == "forward") {
    params.strandprefix=".Forward"
} else {
    params.strandprefix=".Reverse"
}

process Coverage {

	tag "Prefix: $coverage_prefix"
	publishDir "${params.output}/Coverage/wigs",mode:'copy'

	input:
	file inferred_strand from inferred_strand_coverage
	set val(coverage_prefix), file(sorted_coverage_bam), file(sorted_bam_index) from coverage_bam_inputs
	file chr_sizes from chr_sizes

	output:
    file("${coverage_prefix}${params.strandprefix}.wig") into wig_files_temp

	shell:
	'''
	export coverage_strand_rule=$(cat !{inferred_strand})
	if [ $coverage_strand_rule == "none" ] ; then
		!{params.python} $(which !{params.bam2wig}) -s !{chr_sizes} -i !{sorted_coverage_bam} -t !{params.bam2wig_depth_thres} -o !{coverage_prefix}
	else
		!{params.python} $(which !{params.bam2wig}) -s !{chr_sizes} -i !{sorted_coverage_bam} -t !{params.bam2wig_depth_thres} -o !{coverage_prefix} -d "${coverage_strand_rule}"
	fi
	'''
}

/*
 * Step 5b: Wig to BigWig
 */

wig_files_temp
  .flatten()
  .map{ file -> tuple(get_prefix(file), file) }
  .set{ wig_files }

process WigToBigWig {

	tag "Prefix: $wig_prefix"
	publishDir "${params.output}/Coverage/BigWigs",mode:'copy'

	input:
	set val(wig_prefix), file(wig_file) from wig_files
	file chr_sizes from chr_sizes

	output:
	file "*.bw" into coverage_bigwigs

	shell:
	'''
	!{params.wigToBigWig} -clip !{wig_file} !{chr_sizes} !{wig_prefix}.bw
	'''

}

coverage_bigwigs
    .map { f -> tuple(get_read_type(f), f) }
    .groupTuple()
	.into{ mean_coverage_bigwigs;full_coverage_bigwigs }

/*
 * Step 5c: Mean Coverage
 */

process MeanCoverage {
	
	tag "Strand: ${read_type}"
	publishDir "${params.output}/Coverage/mean",'mode':'copy'

	input:
	file inferred_strand_file from inferred_strand_mean_coverage
	set val(read_type), file(mean_coverage_bigwig) from mean_coverage_bigwigs
	file chr_sizes from chr_sizes

	output:
	file "mean*.bw" into mean_bigwigs, expressed_regions_mean_bigwigs

	shell:
	'''
	export coverage_strand_rule=$(cat !{inferred_strand_file})
	if [ $coverage_strand_rule == "none" ] ; then
		!{params.wiggletools} write mean.wig mean !{mean_coverage_bigwig}
		!{params.wigToBigWig} mean.wig !{chr_sizes} mean.bw
	elif [ !{read_type} == "Forward" ]; then
		!{params.wiggletools} write mean.forward.wig mean !{mean_coverage_bigwig}
		!{params.wigToBigWig} mean.forward.wig !{chr_sizes} mean.forward.bw
    else
		!{params.wiggletools} write mean.reverse.wig mean !{mean_coverage_bigwig}
		!{params.wigToBigWig} mean.reverse.wig !{chr_sizes} mean.reverse.bw
	fi
	'''
}

//for hg38, hg19, and mm10, step 6 is enabled by params.step6 = true during Define Reference Paths/Scripts + Reference Dependent Parameters
if (params.step6) {

	/*
	 * Step 6: txQuant
	 */

  process TXQuant {

		tag "Prefix: $salmon_input_prefix"
		publishDir "${params.output}/Salmon_tx/${salmon_input_prefix}",mode:'copy'

		input:
		file salmon_index from salmon_index
		set val(salmon_input_prefix), file(salmon_inputs) from salmon_inputs

		output:
		file "${salmon_input_prefix}/*"
		file("${salmon_input_prefix}_quant.sf") into salmon_quants

		shell:
		salmon_index_prefix = "${params.salmon_prefix}"
		if (params.sample == "single") {
			sample_command = "-r ${salmon_input_prefix}.f*q*"
			if (params.strand == "unstranded" ) {
				salmon_strand = "U"
			} else if (params.strand == "forward" ) {
				salmon_strand = "SF"
			} else { // reverse
				salmon_strand = "SR"
			}
		} else { // paired
			sample_command = "-1 ${salmon_input_prefix}_1.f*q* -2 ${salmon_input_prefix}_2.f*q*"
			if (params.strand == "unstranded" ) {
				salmon_strand = "IU"
			} else if (params.strand == "forward" ) {
				salmon_strand = "ISF"
			} else { // reverse
				salmon_strand = "ISR"
			}
		}
		'''
		!{params.salmon} quant \
		-i !{salmon_index} \
		-p !{task.cpus} \
		-l !{salmon_strand} \
		!{sample_command} \
		-o !{salmon_input_prefix}
		cp !{salmon_input_prefix}/quant.sf !{salmon_input_prefix}_quant.sf
		'''
	}
}

/*
 * Construct the Counts Objects Input Channel
 */

count_objects_bam_files // this puts sorted.bams and bais into the channel
  .flatten()
  .mix(count_objects_quality_reports) //this puts sample_XX_summary.txt files into the channel
  .mix(count_objects_quality_metrics) // this puts sample_XX_fastqc_data.txt into the channel
  .mix(alignment_summaries) // this puts sample_XX_align_summary.txt into the channel
  .mix(create_counts_gtf) // this puts gencode.v25.annotation.gtf file into the channel
  .mix(sample_counts) // !! this one puts sample_05_Gencode.v25.hg38_Exons.counts and sample_05_Gencode.v25.hg38_Genes.counts into the channel
  .mix(regtools_outputs) // !! this one includes the missing *_junctions_primaryOnly_regtools.count files for the CountObjects process
  .flatten()
  .toSortedList()
  .set{ counts_objects_channel }

if (params.reference_type == "human" || params.reference_type == "mouse") {

	counts_objects_channel
	  .mix(salmon_quants)
	  .flatten()
	  .toSortedList()
	  .set{counts_objects_channel_1}
} else {

	counts_objects_channel
	  .set{counts_objects_channel_1}
}
if (params.ercc) {
		
	counts_objects_channel_1
	  .mix(ercc_abundances)
	  .flatten()
	  .toSortedList()
	  .set{ counts_inputs }
} else {

	counts_objects_channel_1
	  .flatten()
	  .toSortedList()
	  .set{ counts_inputs }
}

/*
 * Construct the Annotation Input Channel
 */

//  Mix with reference-dependent annotation info
if(params.reference == "hg19") {
  junction_annotation_ensembl
    .mix(junction_annotation_genes)
    .mix(junction_annotation_gencode)
    .mix(feature_to_tx_gencode)
    .mix(feature_to_tx_ensembl)
    .toSortedList()
    .set{counts_annotations}
} else if (params.reference == "hg38") {
  junction_annotation_ensembl
    .mix(junction_annotation_genes)
    .mix(junction_annotation_gencode)
    .mix(feature_to_tx_gencode)
    .mix(feature_to_tx_ensembl)
    .mix(exon_maps_by_coord_hg38)
    .toSortedList()
    .set{counts_annotations}
} else if (params.reference == "mm10") {
  junction_annotation_ensembl
    .mix(junction_annotation_gencode)
    .toSortedList()
    .set{counts_annotations}
} else {  // rat
  junction_annotation_ensembl
    .toSortedList()
    .set{counts_annotations}
}

/*
 * Step 7a: Create Count Objects
*/

process CountObjects {

	publishDir "${params.output}/Count_Objects",'mode':'copy'

	input:
	file counts_input from counts_inputs
	file counts_annotation from counts_annotations
	file create_counts from create_counts
	file ercc_actual_conc from ercc_actual_conc
	file counts_sample_manifest from counts_samples_manifest

	output:
	file "*"

	shell:
	if (params.ercc) {
		ercc_bool = "TRUE"
	} else {
		ercc_bool = "FALSE"
	}
	if (params.sample == "paired") {
		counts_pe = "TRUE"
	} else {
		counts_pe = "FALSE"
	}
	if (params.strand == "unstranded") {
		counts_strand = "-s FALSE"
    } else {
        counts_strand = "-s " + params.strand
	}
	'''
  !{params.Rscript} !{create_counts} -o !{params.reference} -m ./ -e !{params.experiment} -p !{params.prefix} -l !{counts_pe} -c !{ercc_bool} -t !{task.cpus} !{counts_strand}
	'''
}

if (params.fullCov) {

	full_coverage_bams
	  .flatten()
	  .mix(full_coverage_bigwigs)
	  .flatten()
	  .toSortedList()
	  .set{ full_coverage_inputs }
  
	/*
	 * Step 7b: Create Full Coverage Objects
	 */

    process CoverageObjects {

        publishDir "${params.output}/Coverage_Objects",'mode':'copy'

        input:
            file fullCov_script from fullCov_script
            file fullCov_samples_manifest from fullCov_samples_manifest
            file full_coverage_input from full_coverage_inputs
            file inferred_strand_R_object from inferred_strand_objects

        output:
            file "*"
    
        shell:
            if (params.sample == "paired") {
              coverage_pe = "TRUE"
            } else {
              coverage_pe = "FALSE"
            }
            '''
            !{params.Rscript} !{fullCov_script} -o !{params.reference} -m . -e !{params.experiment} -p !{params.prefix} -l !{coverage_pe} -f TRUE -c !{task.cpus}
            '''
    }
}

if (params.step8) {

	/*
	 * Step 8: Call Variants
	*/

	variant_calls_bam
	  .combine(snvbed)
	  .combine(variant_assembly)
	  .set{ variant_calls }

	process VariantCalls {
 
		tag "Prefix: $variant_bams_prefix"
		publishDir "${params.output}/Variant_Calls",'mode':'copy'

		input:
		set val(variant_bams_prefix), file(variant_calls_bam_file), file(variant_calls_bai), file(snv_bed), file(variant_assembly_file) from variant_calls

		output:
		file "${variant_bams_prefix}.vcf.gz" into compressed_variant_calls
		file "${variant_bams_prefix}.vcf.gz.tbi" into compressed_variant_calls_tbi

		shell:
		'''
		!{params.samtools} mpileup -l !{snv_bed} -AB -q !{params.samtools_min_map_q} -Q !{params.samtools_min_base_q} -d !{params.samtools_max_depth} -uf !{variant_assembly_file} !{variant_calls_bam_file} -o !{variant_bams_prefix}_tmp.vcf
		!{params.bcftools} call -mv -Oz !{variant_bams_prefix}_tmp.vcf > !{variant_bams_prefix}.vcf.gz
		!{params.tabix} -p vcf !{variant_bams_prefix}.vcf.gz
		'''
	}

	compressed_variant_calls
	  .collect()
	  .set{ collected_variant_calls }

	compressed_variant_calls_tbi
	  .collect()
	  .set{ collected_variant_calls_tbi }


	/*
	 * Step 8b: Merge Variant Calls
	 */

	process VariantsMerge {
	
		tag "Multi-sample vcf creation"
		publishDir "${params.output}/Merged_Variants",'mode':'copy'

		input:
		file collected_variants from collected_variant_calls
		file collected_variants_tbi from collected_variant_calls_tbi

		output:
		file "*"

		shell:
		'''
		!{params.bcftools} merge !{collected_variants} | !{params.bgzip} -c > mergedVariants.vcf.gz
		'''
	}
}

/*
 * Step 9: Expressed Regions
 */

process ExpressedRegions {
 
  tag "Sample: $expressed_regions_mean_bigwig"
  publishDir "${params.output}/Expressed_Regions",mode:'copy'

  input:
    file expressed_regions_script from expressed_regions_script
    file chr_sizes from chr_sizes
    file expressed_regions_mean_bigwig from expressed_regions_mean_bigwigs

  output:
    file "*" 

  shell:
    '''
    for meanfile in ./mean*.bw
    do
      !{params.Rscript} !{expressed_regions_script} \
      -m ${meanfile} \
      -o . \
      -i !{chr_sizes} \
      -c !{task.cpus}
    done
    '''
}
