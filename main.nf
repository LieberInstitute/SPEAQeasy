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
	-   II: Download GENCODE GTF File
	- IIIa: Download Transcript Fasta
	- IIIb: Build Salmon or Kallisto Index
	Sample Processing:
	-   A: Input preprocessing
	-   B: ERCC Quality Analysis (Optional)
	-   1: FastQC Quality Analysis
	-  2a: File Trimming (Sample-dependent)
	-  2b: FastQC on trimmed samples
	-  3a: Hisat2 Alignment
	-  3b: Convert Sam to Bam
	-  4a: Feature Counts
	-  4b: Primary Alignments
	-  4c: Junctions
	-  5a: Coverage
	-  5b: WigtoBigWig
	-  5c: Mean Coverage
	-   6: Transcript Quantification (Salmon or Kallisto)
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
      unstranded <- uses pipeline inferencing
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
	--input			The directory containing samples.manifest, the file describing the input FASTQ files. Defaults to "./input" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--output		Defines the output folder for the files. Defaults to "./results" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--unalign		Give the option to not algin the reads against a reference in HISAT step. Defaults to false 
	-----------------------------------------------------------------------------------------------------------------------------------
	--annotation	Path to the folder containing pipeline annotations. Defaults to "./Annotations" (relative to the repository)
	-----------------------------------------------------------------------------------------------------------------------------------
	--ercc			Flag to enable ERCC quantification with Kallisto
	-----------------------------------------------------------------------------------------------------------------------------------
	--fullCov		Flag to perform full coverage in step 7b
	-----------------------------------------------------------------------------------------------------------------------------------
	--small_test	Runs the pipeline as a small test run on sample files located in the test folder
	-----------------------------------------------------------------------------------------------------------------------------------
    --force_trim    Include to perform trimming on all inputs, rather than just those failing QC due to detected adapter content
    -----------------------------------------------------------------------------------------------------------------------------------
    --use_salmon  Include this flag to perform transcript quantification with salmon (which output objects will reflect), rather than just kallisto
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
params.output = "${workflow.projectDir}/results"
params.scripts = "${workflow.projectDir}/scripts"
params.ercc = false
params.fullCov = false
params.small_test = false
params.force_trim = false
params.use_salmon = false
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
// Define Reference Paths/Scripts + Reference-dependent Parameters
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Variant calling is only enabled for human
if (params.reference_type == "human") {
    params.step8 = true
    snvbed = Channel.fromPath("${params.annotation}/Genotyping/common_missense_SNVs_${params.reference}.bed")
} else {
    params.step8 = false
}

if (params.reference == "hg38") {
  params.anno_version = params.gencode_version_human
  params.anno_suffix = params.reference + '_gencode_v' + params.anno_version + '_main'
	
	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh38.primary_assembly.genome.fa.gz"

	// Step 4: gencode gtf
	params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.annotation.gtf.gz"

	// Step 6: transcript quantification
	params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.transcripts.fa.gz"

    // Our extra exon annotation if user defaults to gencode release 25
    exon_maps_by_coord_hg38 = Channel.fromPath("${params.annotation}/junction_txdb/exonMaps_by_coord_hg38_gencode_v25.rda")

} else if (params.reference == "hg19") {
  params.anno_version = params.gencode_version_human
  params.anno_suffix = params.reference + '_gencode_v' + params.anno_version + 'lift37_main'
	
	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
	
	// Step 4: gencode gtf
	params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/gencode.v${params.anno_version}lift37.annotation.gtf.gz"

	// Step 6: transcript quantification
	params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/gencode.v${params.anno_version}lift37.transcripts.fa.gz"

} else if (params.reference == "mm10") {
  params.anno_version = params.gencode_version_mouse
  params.anno_suffix = params.reference + '_gencode_' + params.anno_version + '_main'

	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/GRCm38.primary_assembly.genome.fa.gz"

	// Step 4: gencode gtf
	params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.annotation.gtf.gz"

	// Step 6: transcript quantification
	params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.transcripts.fa.gz"

}
if (params.reference == "rn6") {
    params.anno_version = params.ensembl_version_rat
    params.anno_suffix = params.reference + '_ensembl_' + params.anno_version + '_main'
    
    // Step 3: hisat2
    params.fa_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"

    // Step 4: ensembl gtf
    params.gtf_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.${params.anno_version}.gtf.gz"

    // Step 6: transcript quantification
    params.tx_fa_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Utilities for retrieving info from filenames
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_prefix(f) {
  //  Remove these regardless of position in the string
  String blackListAny = "(_1|_2)*_summary|(_1|_2)*_fastqc_data|_trimmed|_untrimmed|_unpaired|_paired|_hisat_out"
  
  f.name.toString()
   .replaceAll("_1\\.", ".")
   .replaceAll("_2\\.", ".")
   .replaceAll("\\.Forward", "_Forward")
   .replaceAll("\\.Reverse", "_Reverse")
   .tokenize('.')[0]
   .replaceAll(blackListAny, "")
}

def get_read_type(f) {
  baseName = f.name.toString().tokenize('.')[0]
  if (baseName.endsWith('_Forward')) {
    return('forward')
  } else if (baseName.endsWith('_Reverse')) {
    return('reverse')
  } else {
    return('unstranded')
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
summary['Input']			   = params.input
summary['scripts']      = params.scripts
summary['Experiment'] = params.experiment
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


// Uses "storeDir" to download assembly file if required, and otherwise output the cached file
// if it already exists
process pullGENCODEassemblyfa {

    tag "Downloading Assembly FA File: ${baseName}"
    storeDir "${params.annotation}/reference/${params.reference}/assembly/fa"
        
    output:
        file "${out_fasta}" into reference_fasta, variant_assembly, annotation_assembly

    shell:
        //  Name of the primary assembly fasta
        baseName = file("${params.fa_link}").getName() - ".gz"
        
        //  Which fasta to use for this pipeline execution instance
        if (params.gencode_build == "main") {
            out_fasta = "assembly_${params.anno_suffix}.fa"
        } else {
            out_fasta = baseName
        }
        '''
        #  Pull and unzip primary assembly fasta, if this hasn't been done. Otherwise
        #  symbolically link the fasta into the working directory
        if [ -f !{params.annotation}/reference/!{params.reference}/assembly/fa/!{baseName} ]; then
            ln -s !{params.annotation}/reference/!{params.reference}/assembly/fa/!{baseName} !{baseName}
        else
            wget "!{params.fa_link}"
            gunzip "!{baseName}.gz"
        fi
        
        #######################################################################
        #  Create the "main" fasta of canonical seqs only
        #######################################################################
        
        #  Determine how many chromosomes/seqs to keep
        if [ !{params.reference_type} == "human" ]; then
            num_chrs=25
        elif [ !{params.reference} == "mm10" ]; then
            num_chrs=22
        else
            num_chrs=23
        fi
        
        #  Find the line of the header for the first extra contig (to not
        #  include in the "main" annotation fasta
        first_bad_line=$(grep -n ">" !{baseName} | cut -d : -f 1 | paste -s | cut -f $(($num_chrs + 1)))
        
        #  Make a new file out of all the lines up and not including that
        sed -n "1,$(($first_bad_line - 1))p;${first_bad_line}q" !{baseName} > assembly_!{params.anno_suffix}.fa
        '''
}

/*
 * Step Ib: Build HISAT Index
 */

// Uses "storeDir" to build HISAT2 index only when it doesn't exist, and output the cached
// files if they do already exist
process buildHISATindex {
		
  tag "Building HISAT2 Index: ${hisat_prefix}"
  storeDir "${params.annotation}/reference/${params.reference}/assembly/index"

  input:
    file reference_fasta

  output:
    file("${hisat_prefix}.*") into hisat_index

  shell:
    hisat_prefix = "hisat2_assembly_${params.anno_suffix}"
    '''
    !{params.hisat2build} -p !{task.cpus} !{reference_fasta} !{hisat_prefix}
    '''
}


/*
 * Step II: GENCODE GTF Download
 */

// Uses "storeDir" to download gtf only when it doesn't exist, and output the cached
// file if it does already exist
process pullGENCODEgtf {

    tag "Downloading GTF File: ${baseName}"
    storeDir "${params.annotation}/RSeQC/${params.reference}/gtf"

    output:
        file "${out_gtf}" into create_counts_gtf, gencode_feature_gtf, annotation_gtf
        file "${baseName}"

    shell:
        baseName = file("${params.gtf_link}").getName() - ".gz"

        //  Which fasta to use for this pipeline execution instance
        if (params.gencode_build == "main") {
            out_gtf = "transcripts_${params.anno_suffix}.gtf"
        } else {
            out_gtf = baseName
        }
        '''
        #  Pull and unzip transcript gtf
        wget "!{params.gtf_link}"
        gunzip "!{baseName}.gz"

        #  Create the "main" transcripts gtf excluding non-canonical contigs
        grep -E "^(chr|#|[[:digit:]]).*" !{baseName} > transcripts_!{params.anno_suffix}.gtf
        '''
}

process BuildAnnotationObjects {

  storeDir "${params.annotation}/junction_txdb"
  
  input:
      file annotation_assembly
      file annotation_gtf
      file build_ann_script from file("${params.scripts}/build_annotation_objects.R")
      
  output:
      file "junction_annotation_${params.anno_suffix}.rda" into junction_annotation
      file "feature_to_Tx_${params.anno_suffix}.rda" optional true into feature_to_tx_gencode
      file "chrom_sizes_${params.anno_suffix}" into chr_sizes
      
  shell:
      '''
      !{params.Rscript} !{build_ann_script} -r !{params.reference} -v !{params.anno_version} -t "main"
      '''
}

/*
 * Step IIIa: GENCODE TX FA Download
 */
		
// Uses "storeDir" to download files only when they don't exist, and output the cached
// files if they do already exist
process pullGENCODEtranscripts {
			
    tag "Downloading TX FA File: ${baseName}"
    storeDir "${params.annotation}/reference/${params.reference}/transcripts/fa"

    output:
        file baseName into transcript_fa

    shell:
        baseName = file("${params.tx_fa_link}").getName() - ".gz"
        '''
        wget !{params.tx_fa_link}
        gunzip !{baseName}.gz
        '''
}

/*
 * Step IIIb: Build transcript index for Salmon or Kallisto
 */

// Uses "storeDir" to build the index only if the built file is not present; outputs
// this cached file otherwise
if (params.use_salmon) {
    process buildSALMONindex {
    
        tag "Building Salmon Index: salmon_index_${params.anno_suffix}"
        storeDir "${params.annotation}/reference/${params.reference}/transcripts/salmon"
    
        input:
          file transcript_fa
    
        output:
          file("salmon_index_${params.anno_suffix}") into salmon_index
          file("build_salmon_index.log")
    
        script:
          if (params.reference == "rn6") {
              gencode_flag = ""
          } else {
              gencode_flag = "--gencode"
          } 
          """
          ${params.salmon} index \
              -t $transcript_fa \
              -i salmon_index_${params.anno_suffix} \
              -p $task.cpus \
              $gencode_flag \
              -k ${params.salmon_min_read_len}
          cp .command.log build_salmon_index_${params.anno_suffix}.log
          """
    }
}

process BuildKallistoIndex {

    storeDir "${params.annotation}/reference/${params.reference}/transcripts/kallisto"
    
    input:
        file transcript_fa

    output:
        file "kallisto_index_${params.anno_suffix}" into kallisto_index
        file "build_kallisto_index_${params.anno_suffix}.log"

    shell:
        '''
        !{params.kallisto} index -i kallisto_index_!{params.anno_suffix} !{transcript_fa}
        cp .command.log build_kallisto_index_!{params.anno_suffix}.log
        '''
}


//  Step A: Merge files as necessary, rename files based on sample ids provided
//          in the manifest, and create a new manifest for internal use based
//          on these changes

process PreprocessInputs {
  
  publishDir "${params.output}", mode:'copy', pattern:'*.log'

  input:
    file original_manifest from file("${params.input}/samples.manifest")
    file merge_script from file("${params.scripts}/preprocess_inputs.R")

  output:
    file "*.f*q*" into merged_inputs_flat
    file "samples_processed.manifest" into strandness_manifest
    file "preprocess_inputs.log"
  
  shell:
    '''
    !{params.Rscript} !{merge_script}
    cp .command.log preprocess_inputs.log
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
        .flatten()
        .map{file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
        .set{ temp_inputs }
}

//  Copy the processed input channel to channels used by dependent processes
if (params.ercc) { 
  temp_inputs.into{ strandness_inputs; fastqc_untrimmed_inputs; trimming_fqs; tx_quant_inputs; ercc_inputs }
} else {
  temp_inputs.into{ strandness_inputs; fastqc_untrimmed_inputs; trimming_fqs; tx_quant_inputs }
}

process InferStrandness {

    tag "${prefix}"
    publishDir "${params.output}/InferStrandness", mode:'copy'
    
    input:
        file infer_strand_R from file("${params.scripts}/infer_strand.R")
        file infer_strand_sh from file("${params.scripts}/infer_strand.sh")
        file kallisto_index
        set val(prefix), file(fq_file) from strandness_inputs
        
    output:
        file "${prefix}_infer_strand.log"
        file "${prefix}_strandness_pattern.txt" into strandness_patterns
        
    shell:
        '''
        bash !{infer_strand_sh} \
            !{params.sample} \
            !{params.num_reads_infer_strand} \
            !{params.kallisto} \
            !{task.cpus} \
            !{params.kallisto_len_mean} \
            !{params.kallisto_len_sd} \
            !{params.strand} \
            !{params.Rscript}
        
        cp .command.log !{prefix}_infer_strand.log
        '''
}

process CompleteManifest {

    publishDir "${params.output}", mode:'copy'
    
    input:
        file strandness_files from strandness_patterns.collect()
        file strandness_manifest
        file manifest_script from file("${params.scripts}/complete_manifest.R")
        
    output:
        file "samples_complete.manifest" into complete_manifest_ercc, complete_manifest_feature, complete_manifest_junctions, complete_manifest_cov, complete_manifest_counts, complete_manifest_hisat, complete_manifest_quant, complete_manifest_fullcov
        
    shell:
        '''
        !{params.Rscript} !{manifest_script}
        '''
}
    
/*
 * Step B: Run the ERCC process if the --ercc flag is specified
 */
 
if (params.ercc) {

  process ERCC {
		
    tag "Prefix: $prefix"
    publishDir "${params.output}/ercc/${prefix}",'mode':'copy'

    input:
      file erccidx from file("${params.annotation}/ERCC/ERCC92.idx")
      set val(prefix), file(ercc_input) from ercc_inputs
      file complete_manifest_ercc

    output:
      file "${prefix}_abundance.tsv" into ercc_abundances
      file "ercc_${prefix}.log"

    shell:
      if (params.sample == "single") {
          kallisto_flags = "--single -l ${params.kallisto_len_mean} -s ${params.kallisto_len_sd}"
      } else {
          kallisto_flags = ""
      }
      '''
      #  Find this sample's strandness and determine kallisto flags to set
      strand=$(cat samples_complete.manifest | grep !{prefix} | awk -F ' ' '{print $NF}')
      if [ $strand == 'forward' ]; then
          kallisto_strand=" --fr-stranded"
      elif [ $strand == 'reverse' ]
          kallisto_strand=" --rf-stranded"
      else
          kallisto_strand=""
      fi
      
      !{params.kallisto} quant -i !{erccidx} -t !{task.cpus} !{kallisto_flags} -o . $kallisto_strand !{ercc_input}
      cp abundance.tsv !{prefix}_abundance.tsv
      cp .command.log ercc_!{prefix}.log
      '''
  }
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
	${params.fastqc} -t $task.cpus $fastqc_untrimmed_input --extract
	$copy_command
	$data_command
	"""
}

//  Combine FASTQ files and FastQC result summaries for each sample, to form the input channel for Trimming
if (params.sample == "single") {

    quality_reports
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .join(trimming_fqs)
        .ifEmpty{ error "All files (fastQC summaries on untrimmed inputs, and the FASTQs themselves) missing from input to trimming channel." }
        .set{ trimming_inputs }
        
} else { // paired

    quality_reports
        .flatten()
        .map{ file -> tuple(get_prefix(file), file) }
        .groupTuple()
        .join(trimming_fqs)
        .ifEmpty{ error "All files (fastQC summaries on untrimmed inputs, and the FASTQs themselves) missing from input to trimming channel." }
        .set{ trimming_inputs }
}    

/*
 * Step 2a: Trim samples (dependent on user-chosen settings)
 */
 
process Trimming {

    tag "Prefix: $fq_prefix"
    publishDir "${params.output}/Trimming",mode:'copy'

    input:
        set val(fq_prefix), file(fq_summary), file(fq_file) from trimming_inputs

    output:
        file "${fq_prefix}_trimmed*.fastq" optional true into trimmed_fastqc_inputs
        file "${fq_prefix}*.f*q*" into trimming_outputs

    shell:
        file_ext = get_file_ext(fq_file[0])
        if (params.sample == "single") {
            output_option = "${fq_prefix}_trimmed.fastq"
            trim_mode = "SE"
            adapter_fa_temp = params.adapter_fasta_single
            trim_clip = params.trim_clip_single
        } else {
            output_option = "${fq_prefix}_trimmed_paired_1.fastq ${fq_prefix}_unpaired_1.fastq ${fq_prefix}_trimmed_paired_2.fastq ${fq_prefix}_unpaired_2.fastq"
            trim_mode = "PE"
            adapter_fa_temp = params.adapter_fasta_paired
            trim_clip = params.trim_clip_paired
        }
        '''
        #  Determine whether to trim the FASTQ file(s). This is done if the user
        #  adds the --force_trim flag, or if fastQC adapter content metrics fail
        #  for at least one FASTQ
        if [ "!{params.force_trim}" == true ]; then
            do_trim=true
        elif [ "!{params.sample}" == "single" ]; then
            if [ $(grep "Adapter Content" !{fq_summary} | cut -f 1)  == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        else
            result1=$(grep "Adapter Content" !{fq_prefix}_1_summary.txt | cut -c1-4)
            result2=$(grep "Adapter Content" !{fq_prefix}_2_summary.txt | cut -c1-4)
            if [ $result1 == "FAIL" ] || [ $result2 == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        fi
        
        #  Run trimming if required
        if [ "$do_trim" == true ]; then
            #  This solves the problem of trimmomatic and the adapter fasta
            #  needing hard paths, even when on the PATH.
            if [ !{params.use_long_paths} == "true" ]; then
                trim_jar=!{params.trimmomatic}
                adapter_fa=!{adapter_fa_temp}
            else
                trim_jar=$(which !{params.trimmomatic})
                adapter_fa=$(which !{adapter_fa_temp})
            fi

            java -Xmx512M \
                -jar $trim_jar \
                !{trim_mode} \
                -threads !{task.cpus} \
                -phred33 \
                *.f*q* \
                !{output_option} \
                ILLUMINACLIP:$adapter_fa:!{trim_clip} \
                LEADING:!{params.trim_lead} \
                TRAILING:!{params.trim_trail} \
                SLIDINGWINDOW:!{params.trim_slide_window} \
                MINLEN:!{params.trim_min_len}
        else
            #  Otherwise rename files (signal to nextflow to output these files)
            if [ "!{params.sample}" == "single" ]; then
                mv !{fq_prefix}!{file_ext} !{fq_prefix}_untrimmed!{file_ext}
            else
                mv !{fq_prefix}_1!{file_ext} !{fq_prefix}_untrimmed_1!{file_ext}
                mv !{fq_prefix}_2!{file_ext} !{fq_prefix}_untrimmed_2!{file_ext}
            fi
        fi
        
        cp .command.log trimming_!{fq_prefix}.log
        '''
}


/*
 * Step 2b: Run FastQC Quality Check on Trimmed Files
 */

process QualityTrimmed {

	tag "Prefix: $fastq_name"
	publishDir "${params.output}/FastQC/Trimmed",'mode':'copy'

	input:
	file fastqc_trimmed_input from trimmed_fastqc_inputs

	output:
	file "*"

	script:
    fastq_name = get_prefix(fastqc_trimmed_input[0])
	"""
	$params.fastqc -t $task.cpus $fastqc_trimmed_input --extract
	"""
}

/*
 * Step 3a: Hisat Alignment
 */

trimming_outputs
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .ifEmpty{ error "Input channel to HISAT is empty" }
    .groupTuple()
    .set{ hisat_inputs }
        
if (params.sample == "single") {

    process SingleEndHISAT {

        tag "Prefix: $prefix"
        publishDir "${params.output}/HISAT2_out",mode:'copy'

        input:
            file hisat_index
            file complete_manifest_hisat
            set val(prefix), file(single_hisat_input) from hisat_inputs

        output:
            file "*_hisat_out.sam" into hisat_output
            file "*"
            file "*_align_summary.txt" into alignment_summaries

        shell:
            '''
            #  Find this sample's strandness and determine strand flag
            strand=$(cat samples_complete.manifest | grep !{prefix} | awk -F ' ' '{print $NF}')
            if [ ${strand} == "unstranded" ]; then
                hisat_strand=""
            elif [ ${strand} == "forward" ]; then
                hisat_strand="--rna-strandness F"
            else
                hisat_strand="--rna-strandness R"
            fi
            
            #  Run Hisat2
            !{params.hisat2} \
                -p !{task.cpus} \
                -x !{params.annotation}/reference/!{params.reference}/assembly/index/hisat2_assembly_!{params.anno_suffix} \
                -U !{single_hisat_input} \
                -S !{prefix}_hisat_out.sam \
                ${hisat_strand} \
                --phred33 \
                --min-intronlen !{params.min_intron_len} \
                2> !{prefix}_align_summary.txt
            '''
    }
} else { // sample is paired-end

    process PairedEndHISAT {

        tag "Prefix: $prefix"
        publishDir "${params.output}/HISAT2_out",'mode':'copy'

        input:
            file hisat_index
            file complete_manifest_hisat
            set val(prefix), file(fq_files) from hisat_inputs

        output:
            file "*_hisat_out.sam" into hisat_output
            file "*"
            file "*_align_summary.txt" into alignment_summaries

        shell:
            if (params.unalign) {
                unaligned_opt = "--un-conc ${prefix}_discordant.fastq"
            } else {
                unaligned_opt = ""
            }
            '''
            #  Find this sample's strandness and determine strand flag
            strand=$(cat samples_complete.manifest | grep !{prefix} | awk -F ' ' '{print $NF}')
            if [ ${strand} == "unstranded" ]; then
                hisat_strand=""
            elif [ ${strand} == "forward" ]; then
                hisat_strand="--rna-strandness FR"
            else
                hisat_strand="--rna-strandness RF"
            fi
            
            #  If this sample had unpaired trimming outputs, include them
            if [ -f !{prefix}_unpaired_1.fastq ]; then
                unpaired_opt='-U !{prefix}_unpaired_1.fastq,!{prefix}_unpaired_2.fastq'
            else
                unpaired_opt=''
            fi
            
            #  Run Hisat2
            !{params.hisat2} \
                -p !{task.cpus} \
                -x !{params.annotation}/reference/!{params.reference}/assembly/index/hisat2_assembly_!{params.anno_suffix} \
                -1 !{prefix}*trimmed*_1.f*q* \
                -2 !{prefix}*trimmed*_2.f*q* \
                -S !{prefix}_hisat_out.sam \
                ${unpaired_opt} \
                ${hisat_strand} \
                --phred33 \
                --min-intronlen !{params.min_intron_len} \
                !{unaligned_opt} \
                2> !{prefix}_align_summary.txt
            '''
    }
}

hisat_output
    .flatten()
    .map{ file -> tuple(get_prefix(file), file) }
    .set{ sam_to_bam_inputs }

/*
 * Step 3b: Sam to Bam 
 */

process SamtoBam {
	
	tag "Prefix: $sam_to_bam_prefix"
	publishDir "${params.output}/HISAT2_out/sam_to_bam",'mode':'copy'

	input:
	set val(sam_to_bam_prefix), file(sam_to_bam_input) from sam_to_bam_inputs

	output:
	set val("${sam_to_bam_prefix}"), file("${sam_to_bam_prefix}*.sorted.bam"), file("${sam_to_bam_prefix}*.sorted.bam.bai") into feature_bam_inputs, alignment_bam_inputs, coverage_bam_inputs, full_coverage_bams, count_objects_bam_files, variant_calls_bam

	script:
	original_bam = "${sam_to_bam_prefix}_accepted_hits.bam"
	sorted_bam = "${sam_to_bam_prefix}_accepted_hits.sorted"
	"""
	${params.samtools} view -bh -F 4 $sam_to_bam_input > $original_bam
	${params.samtools} sort -T temporary -@ $task.cpus $original_bam -o ${sorted_bam}.bam
	${params.samtools} index ${sorted_bam}.bam
	"""
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
        file complete_manifest_feature

    output:
        file "*"
        file "*.counts*" into sample_counts

    script:
        if (params.sample == "single") {
            sample_option = ""
        } else {
            sample_option = "-p"
        }
        feature_out = "${feature_prefix}_${params.anno_suffix}"
 
        """
        #  Find this sample's strandness and determine strand flag
        strand=\$(cat samples_complete.manifest | grep ${feature_prefix} | awk -F ' ' '{print \$NF}')
        if [ \$strand == 'forward' ]; then
            feature_strand=1
        elif [ \$strand == 'reverse' ]; then
            feature_strand=2
        else
            feature_strand=0
        fi
        
        ${params.featureCounts} \
            -s \$feature_strand \
            $sample_option \
            ${params.feat_counts_gene_opts} \
            -T $task.cpus \
            -a $gencode_gtf_feature \
            -o ${feature_out}_Genes.counts \
            $feature_bam

        ${params.featureCounts} \
            -s \$feature_strand \
            $sample_option \
            ${params.feat_counts_exon_opts} \
            -f \
            -T $task.cpus \
            -a $gencode_gtf_feature \
            -o ${feature_out}_Exons.counts \
            $feature_bam

        cp .command.log feature_counts_${feature_prefix}.log
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

    tag "Prefix: $prefix"
    publishDir "${params.output}/Counts/junction",'mode':'copy'

    input:
        file bed_to_juncs_script from file("${params.scripts}/bed_to_juncs.py")
        set val(prefix), file(alignment_bam), file(alignment_index) from primary_alignments
        file complete_manifest_junctions

    output:
        file "*"
        file("*.count") into regtools_outputs

    shell:
        outjxn = "${prefix}_junctions_primaryOnly_regtools.bed"
        outcount = "${prefix}_junctions_primaryOnly_regtools.count"
        '''
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep !{prefix} | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            strand_integer=2
        elif [ $strand == 'reverse' ]; then
            strand_integer=1
        else
            strand_integer=0
        fi
        
        !{params.regtools} junctions extract -m !{params.min_intron_len} -s ${strand_integer} -o !{outjxn} !{alignment_bam}
        python2.7 !{bed_to_juncs_script} < !{outjxn} > !{outcount}
        '''
}

/*
 * Step 5a: Coverage
 */

process Coverage {

    tag "Prefix: $coverage_prefix"
    publishDir "${params.output}/Coverage/wigs",mode:'copy'

    input:
        file complete_manifest_cov
        set val(coverage_prefix), file(sorted_coverage_bam), file(sorted_bam_index) from coverage_bam_inputs
        file chr_sizes from chr_sizes

    output:
        file "${coverage_prefix}*.wig" into wig_files_temp
        file "bam2wig_${coverage_prefix}.log"

    shell:
        '''
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep !{coverage_prefix} | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            if [ !{params.sample} == "paired" ]; then
                strand_flag="-d 1++,1--,2+-,2-+"
            else
                strand_flag="-d ++,--"
            fi
        elif [ $strand == 'reverse' ]; then
            if [ !{params.sample} == "paired" ]; then
                strand_flag="-d 1+-,1-+,2++,2--"
            else
                strand_flag="-d +-,-+"
            fi
        else
            strand_flag=""
        fi

        python2.7 $(which !{params.bam2wig}) \
            -s !{chr_sizes} \
            -i !{sorted_coverage_bam} \
            -t !{params.bam2wig_depth_thres} \
            -o !{coverage_prefix} \
            $strand_flag

        cp .command.log bam2wig_!{coverage_prefix}.log
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
        set val(read_type), file(mean_coverage_bigwig) from mean_coverage_bigwigs
        file chr_sizes from chr_sizes

    output:
        file "mean*.bw" into mean_bigwigs, expressed_regions_mean_bigwigs

    shell:
        '''
        if [ !{read_type} == "unstranded" ] ; then
            !{params.wiggletools} write mean.wig mean !{mean_coverage_bigwig}
            !{params.wigToBigWig} mean.wig !{chr_sizes} mean.bw
        elif [ !{read_type} == "forward" ]; then
            !{params.wiggletools} write mean.forward.wig mean !{mean_coverage_bigwig}
            !{params.wigToBigWig} mean.forward.wig !{chr_sizes} mean.forward.bw
        else
            !{params.wiggletools} write mean.reverse.wig mean !{mean_coverage_bigwig}
            !{params.wigToBigWig} mean.reverse.wig !{chr_sizes} mean.reverse.bw
        fi
        '''
}


/*
 * Step 6: Transcript quantification
 */

if (params.use_salmon) {
    process TXQuantSalmon {
    
        tag "Prefix: $prefix"
        publishDir "${params.output}/Salmon_tx/${prefix}",mode:'copy'
    
        input:
            file salmon_index
            file complete_manifest_quant
            set val(prefix), file(salmon_inputs) from tx_quant_inputs
    
        output:
            file "${prefix}/*"
            file "${prefix}_quant.sf" into tx_quants
            file "tx_quant_${prefix}.log"
    
        shell:
            '''
            #  Find this sample's strandness and determine strand flag
            strand=$(cat samples_complete.manifest | grep !{prefix} | awk -F ' ' '{print $NF}')
            if [ $strand == 'forward' ]; then
                strand_flag="SF"
            elif [ $strand == 'reverse' ]; then
                strand_flag="SR"
            else
                strand_flag="U"
            fi
            
            if [ !{params.sample} == "paired" ]; then
                strand_flag="I$strand_flag"
                sample_flag="-1 !{prefix}_1.f*q* -2 !{prefix}_2.f*q*"
            else
                sample_flag="-r !{prefix}.f*q*"
            fi
            
            !{params.salmon} quant \
                -i !{salmon_index} \
                -p !{task.cpus} \
                -l ${strand_flag} \
                ${sample_flag} \
                -o !{prefix}
    
            cp !{prefix}/quant.sf !{prefix}_quant.sf
            cp .command.log tx_quant_!{prefix}.log
            '''
    }
} else {
    process TXQuantKallisto {
    
        tag "Prefix: $prefix"
        publishDir "${params.output}/kallisto_tx/${prefix}", mode:'copy'
    
        input:
            file kallisto_index
            file complete_manifest_quant
            set val(prefix), file(fastqs) from tx_quant_inputs
    
        output:
            file "*"
            file "${prefix}_abundance.tsv" into tx_quants
            
        shell:
            '''
            #  Find this sample's strandness
            strand=$(cat samples_complete.manifest | grep !{prefix} | awk -F ' ' '{print $NF}')
            if [ $strand == 'forward' ]; then
                strand_flag="--fr-stranded"
            elif [ $strand == 'reverse' ]; then
                strand_flag="--rf-stranded"
            else
                strand_flag=""
            fi
            
            #  Run the quantification step
            if [ !{params.sample} == "paired" ]; then
                fq1=$(ls *_1.f*q*)
                fq2=$(ls *_2.f*q*)
            
                !{params.kallisto} quant -t !{task.cpus} $strand_flag -i !{kallisto_index} -o . $fq1 $fq2
            else
                fq=$(ls *.f*q*)
                
                !{params.kallisto} quant \
                    -t !{task.cpus} \
                    --single \
                    -l !{params.kallisto_len_mean} \
                    -s !{params.kallisto_len_sd} \
                    $strand_flag \
                    -i !{kallisto_index} \
                    -o . $fq
            fi
            
            mv abundance.tsv !{prefix}_abundance.tsv
            cp .command.log tx_quant_!{prefix}.log
            '''
    }
}

/*
 * Construct the Counts Objects Input Channel
 */

count_objects_bam_files // this puts sorted.bams and bais into the channel
  .flatten()
  .mix(count_objects_quality_reports) // this puts sample_XX_summary.txt files into the channel
  .mix(count_objects_quality_metrics) // this puts sample_XX_fastqc_data.txt into the channel
  .mix(alignment_summaries) // this puts sample_XX_align_summary.txt into the channel
  .mix(create_counts_gtf) // this puts gencode.v25.annotation.gtf file into the channel
  .mix(sample_counts) // this puts sample_XX_*_Exons.counts and sample_XX_*_Genes.counts into the channel
  .mix(regtools_outputs) // this puts the *_junctions_primaryOnly_regtools.count files into the channel
  .mix(tx_quants)
  .set{counts_objects_channel_1}

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
if(params.reference == "hg19" || (params.reference == "hg38" && params.anno_version != 25)) {
    junction_annotation
        .mix(feature_to_tx_gencode)
        .toList()
        .set{counts_annotations}
} else if (params.reference == "hg38" && params.anno_version == 25) {
    junction_annotation
        .mix(feature_to_tx_gencode)
        .mix(exon_maps_by_coord_hg38)
        .toList()
        .set{counts_annotations}
} else { // mouse or rat
    junction_annotation
        .set{counts_annotations}
}

/*
 * Step 7a: Create Count Objects
*/

process CountObjects {

    publishDir "${params.output}/Count_Objects",'mode':'copy'

    input:
        file counts_inputs
        file counts_annotations
        file create_counts from file("${params.scripts}/create_count_objects-${params.reference_type}.R")
        file ercc_actual_conc from file("${params.annotation}/ERCC/ercc_actual_conc.txt")
        file complete_manifest_counts

    output:
        file "*"

    shell:
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
        !{params.Rscript} !{create_counts} \
            -o !{params.reference} \
            -e !{params.experiment} \
            -p !{params.prefix} \
            -l !{counts_pe} \
            -c !{params.ercc} \
            -t !{task.cpus} \
            !{counts_strand} \
            -n !{params.use_salmon}

        cp .command.log counts.log
        '''
}

if (params.fullCov) {

	full_coverage_bams
	  .flatten()
	  .mix(full_coverage_bigwigs)
	  .flatten()
	  .collect()
	  .set{ full_coverage_inputs }
  
	/*
	 * Step 7b: Create Full Coverage Objects
	 */

    process CoverageObjects {

        publishDir "${params.output}/Coverage_Objects",'mode':'copy'

        input:
            file fullCov_script from file("${params.scripts}/create_fullCov_object.R")
            file complete_manifest_fullcov
            file full_coverage_input from full_coverage_inputs

        output:
            file "*"
    
        shell:
            if (params.sample == "paired") {
              coverage_pe = "TRUE"
            } else {
              coverage_pe = "FALSE"
            }
            '''
            !{params.Rscript} !{fullCov_script} \
                -o !{params.reference} \
                -e !{params.experiment} \
                -l !{coverage_pe} \
                -c !{task.cpus}
            
            cp .command.log coverage_objects.log
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


	/*
	 * Step 8b: Merge Variant Calls
	 */

	process VariantsMerge {
	
		tag "Multi-sample vcf creation"
		publishDir "${params.output}/Merged_Variants",'mode':'copy'

		input:
		file collected_variants from compressed_variant_calls.collect()
		file collected_variants_tbi from compressed_variant_calls_tbi.collect()

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
    file expressed_regions_script from file("${params.scripts}/step9-find_expressed_regions.R")
    file chr_sizes from chr_sizes
    file expressed_regions_mean_bigwig from expressed_regions_mean_bigwigs

  output:
    file "*" 

  shell:
    // "strand" is used for naming the log file for this execution of the process
    strand = expressed_regions_mean_bigwig.toString().replaceAll("mean.|.bw|bw", "")
    if (strand.length() > 0) {
        strand = '_' + strand
    }
    '''
    for meanfile in ./mean*.bw
    do
      !{params.Rscript} !{expressed_regions_script} \
      -m ${meanfile} \
      -o . \
      -i !{chr_sizes} \
      -c !{task.cpus}
    done
    
    cp .command.log expressed_regions!{strand}.log
    '''
}
