#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
===============================================================================
    SPEAQeasy: an RNA-seq analysis pipeline from LIBD
===============================================================================

RNA-Seq Multi-Input Analysis Pipeline. Nextflow Version: Started December 2017.
 
## Homepage / Documentation
https://github.com/LieberInstitute/SPEAQeasy
 
## Authors
    ## Original Pipeline
    Emily Burke <emily.burke@libd.org>
    Leonardo Collado-Tores <lcolladotor@gmail.com>
    Andrew Jaffe <andrew.jaffe@libd.org>
    BaDoi Phan <badoi.phan@pitt.edu>
 
    ## Nextflow Version
    Nick Eagles <nickeagles77@gmail.com>
    Jacob Leonard <leonard.jacob09@gmail.com>
    Israel Aguilar <iaguilaror@gmail.com>
    Violeta Larios <siedracko@gmail.com>
    
## References
https://github.com/SciLifeLab/NGI-smRNAseq
-------------------------------------------------------------------------------

 Pipeline Overview:

    Preprocessing:
        - Download Assembly FA
	    - Build HISAT Index
        - Download GENCODE GTF File
        - Download Transcript Fasta
        - Build Salmon or Kallisto Index
    Sample Processing:
        - Input preprocessing
        - ERCC Quality Analysis (Optional)
        - FastQC Quality Analysis
        - File Trimming (Sample-dependent)
        - FastQC on trimmed samples
        - Hisat2 Alignment
        - Convert Sam to Bam
        - Feature Counts
        - Primary Alignments
        - Junctions
        - Coverage
        - WigtoBigWig
        - Mean Coverage
        - Transcript Quantification (Salmon or Kallisto)
        - Create Counts Objects
        - Create Coverage Objects
        - Call Variants
        - Merge Called Variants
        - Expressed Regions

===============================================================================
*/

def helpMessage() {
	log.info"""
	=============================================================
    SPEAQeasy: an RNA-seq analysis pipeline from LIBD
	=============================================================
	Usage:
 
	A typical command for running the pipeline is as follows:

	nextflow run main.nf --sample "single" --strand "unstranded" \\
      --reference "hg19" --ercc --fullCov -profile local

	NOTES: The pipeline accepts single or paired end reads. These reads can be
         stranded or unstranded. Base names of FASTQ files can not contain "."
         (e.g. "sample.1.fastq")

  nextflow run main.nf {CORE} {OPTIONS}
  
    {CORE (mandatory)}:
    
    --sample "single"/"paired"
      single <- reads files individually
      paired <- reads files paired
    --reference
      hg38 <- uses human reference hg38
      hg19 <- uses human reference hg19
      mm10 <- uses mouse reference mm10
      rat  <- uses rattus norvegicus reference Rnor_6.0 for annotation release
              < 105, and mRatBN7.2 for annotation release >= 105
      rn6  <- (deprecated) equivalent to "rat"
    --strand "forward"/"reverse"/"unstranded"
      forward    <- asserts reads are forward stranded
      reverse    <- asserts reverse strandness
      unstranded <- asserts no consistent strandness
      
    {OPTIONS}:
    
    --annotation [path] <- the directory to store and check for annotation-
                       related files. Default: "[SPEAQeasy dir]/Annotation".
                       Note only full paths may be used!
    --coverage      <- Include this flag to produce coverage bigWigs and
                       compute expressed genomic regions. These steps are a
                       useful precursor for analyses involving finding
                       differentially expressed regions (DERs). Default: false
    --custom_anno [string] <- Include this flag to state that the directory
                       specified by "--annotation [path]" has user-provided
                       annotation files to use, and objects built from these
                       files should be labelled with the specified string. See
                       the README for specifics. Default: "" (indicating to
                       check for and potentially pull annotation specified by
                       gencode/ensembl version)
    --ercc          <- Flag to enable ERCC quantification with Kallisto
                       (default: false)
    --experiment [string] <- main name for the experiment (default:
                       "Jlab_experiment")
    --fullCov       <- flag to perform full coverage analysis (default: false).
                       Implies the '--coverage' option.
    --help          <- shows this message
    --input [path]  <- the directory containing samples.manifest, the file
                       describing the input FASTQ files. Default:
                       "[SPEAQeasy dir]/input". Note only full paths may be
                       used!
    --keep_unpaired <- include this flag to keep unpaired reads output from
                       trimming paired-end samples, for use in alignment.
                       Default: false, as this can cause issues in downstream
                       tools like FeatureCounts.
    --output [path] <- the directory to place pipeline outputs. Default:
                       "[SPEAQeasy dir]/results". Note only full paths may be used!
    --prefix [string] <- an additional identifier (name) for the experiment
                       (e.g. date, genome)
    --qsva [path]   <- Optional full path to a text file, containing one Ensembl
                       transcript ID per line for each transcript desired in the
                       final transcripts R output object (called `rse_tx`)
    --small_test	  <- use small test files as input, rather than the files
                       specified in samples.manifest in the directory given by
                       "--input [path]". Default: false.
    --strand_mode [mode] <- determines how to handle disagreement between
                       user-asserted and SPEAQeasy-inferred strandness:
                           "accept": warn about disagreement but continue,
                               using SPEAQeasy-inferred strandness downstream
                           "declare": warn about disagreement but continue,
                               using user-asserted strandness downstream
                           "strict": (default) halt with an error if any
                               disagreement occurs (since this often indicates
                               problems in the input data)
    --trim_mode [mode] <- determines the conditions under which trimming occurs:
                          "skip": do not perform trimming on samples
                          "adaptive": [default] perform trimming on samples that
                              have failed the FastQC "Adapter content" metric
                          "force": perform trimming on all samples
    --unalign       <- include this flag to additionally save discordant reads
                       (when using HISAT2) or unmapped reads (when using STAR)
                       to the pipeline outputs (default: false; only concordant
                       .sam files are saved)
    --use_salmon    <- include this flag to perform transcript quantification
                       with salmon instead of the default of kallisto. Default:
                       false
    --use_star      <- include this flag to use STAR during alignment, instead
                       of the default of HISAT2
                       
    The above was a comprehensive list of options specific to SPEAQeasy. As a
    pipeline based on nextflow, SPEAQeasy also accepts any options the
    "nextflow run" command accepts. To print the list of these, type:
    
          nextflow run -h

  ------------------------------------------------------------------------------
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

params.annotation = "${workflow.projectDir}/Annotation"
params.coverage = false
params.custom_anno = ""
params.ercc = false
params.experiment = "Jlab_experiment"
params.fullCov = false
params.keep_unpaired = false
params.output = "${workflow.projectDir}/results"
params.prefix = ""
params.qsva = ""
params.reference = ""
params.sample = ""
params.small_test = false
params.strand = ""
params.strand_mode = "strict"
params.trim_mode = "adaptive"
params.unalign = false
params.use_salmon = false
params.use_star = false

// Re-assign the deprecated rat reference of "rn6" to "rat"
if (params.reference == "rn6") {
    params.reference = "rat"
}

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
if (! (params.reference in ["hg19", "hg38", "mm10", "rat"])) {
    exit 1, "Error: enter 'hg19' or 'hg38' for human, 'mm10' for mouse, or 'rat' for rat as the reference."
}

// Trim mode
if (params.trim_mode != "skip" && params.trim_mode != "adaptive" && params.trim_mode != "force") {
    exit 1, "'--trim_mode' accepts one of three possible arguments: 'skip', 'adaptive', or 'force'."
}

// Strand mode
if (params.strand_mode != "accept" && params.strand_mode != "declare" && params.strand_mode != "strict") {
    exit 1, "'--strand_mode' accepts one of three possible arguments: 'accept', 'declare', or 'strict'."
}

// Keeping unpaired reads that are not produced
if (params.keep_unpaired && params.trim_mode == "skip") {
    exit 1, "You have opted to include unpaired outputs from trimming, but to skip trimming itself. Consider using a different 'trim_mode' or not using the '--keep_unpaired' option."
}

if (params.keep_unpaired && params.use_star) {
    exit 1, "STAR does not support inclusion of unpaired reads. Consider using HISAT2 (default) or remove the '--keep_unpaired' option."
}

// Trying to subset transcripts when they won't be quantified at all
if (params.qsva != "" && params.reference == "rat") {
    println "Warning: ignoring '--qsva' argument, since transcripts are not quantified for rat."
}

// Passing an illegimate file name to '--qsva'
File qsva_file = new File(params.qsva)
if (params.qsva != "" && ! qsva_file.exists()) {
    exit 1, "File passed via '--qsva' argument does not exist."
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
    params.input = "${workflow.projectDir}/test/$params.reference_type/${params.sample}/${params.strand}"
} else {
    params.input = "${workflow.projectDir}/input"
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define reference-dependent variables
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if (params.custom_anno != "") {
    params.anno_suffix = params.custom_anno + "_custom_build"
    params.anno_version = "custom"
} else if (params.reference == "hg38") {
    params.anno_version = params.gencode_version_human
    params.anno_suffix = params.reference + '_gencode_v' + params.anno_version + '_' + params.anno_build

    // Reference assembly fasta, gtf, and transcript fasta
    params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh38.primary_assembly.genome.fa.gz"
    if (params.anno_build == "primary") {
        params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.primary_assembly.annotation.gtf.gz"
    } else {
        params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.annotation.gtf.gz"
    }
    params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/gencode.v${params.anno_version}.transcripts.fa.gz"
} else if (params.reference == "hg19") {
    if (params.anno_build == "primary") {
        print("Warning: use of 'primary' annotation is not supported for hg19, as GENCODE does not provide a primary .gtf file. Continuing with annotation build 'main'.")
    }
    
    params.anno_version = params.gencode_version_human
    params.anno_suffix = params.reference + '_gencode_v' + params.anno_version + 'lift37_main'

    // Reference assembly fasta, gtf, and transcript fasta
    params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
    params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/gencode.v${params.anno_version}lift37.annotation.gtf.gz"
    params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${params.anno_version}/GRCh37_mapping/gencode.v${params.anno_version}lift37.transcripts.fa.gz"

} else if (params.reference == "mm10") {
    params.anno_version = params.gencode_version_mouse
    params.anno_suffix = params.reference + '_gencode_' + params.anno_version + '_' + params.anno_build

    // Reference assembly fasta, gtf, and transcript fasta
    params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/GRCm38.primary_assembly.genome.fa.gz"
    if (params.anno_build == "primary") {
        params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.primary_assembly.annotation.gtf.gz"
    } else {
        params.gtf_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.annotation.gtf.gz"
    }
    params.tx_fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${params.anno_version}/gencode.v${params.anno_version}.transcripts.fa.gz"

} else { // rat
    params.anno_version = params.ensembl_version_rat
    params.anno_suffix = params.reference + '_ensembl_' + params.anno_version + '_' + params.anno_build

    // At Ensembl, the rat genome switches from "Rnor_6.0" to "mRatBN7.2" at and
    // after release 105
    if (Integer.parseInt(params.anno_version) >= 105) {
        rat_genome_name = "mRatBN7.2" 
    } else {
        rat_genome_name = "Rnor_6.0" 
    }
    // Reference assembly fasta, gtf, and transcript fasta
    params.fa_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.${rat_genome_name}.dna.toplevel.fa.gz"
    if (params.anno_build == "primary") {
        params.gtf_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.${rat_genome_name}.${params.anno_version}.gtf.gz"
    } else {
        params.gtf_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/gtf/rattus_norvegicus/Rattus_norvegicus.${rat_genome_name}.${params.anno_version}.chr.gtf.gz"
    }
    params.tx_fa_link = "ftp://ftp.ensembl.org/pub/release-${params.anno_version}/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.${rat_genome_name}.cdna.all.fa.gz"
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    Utilities for retrieving info from filenames
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_ids(row) {
    if (params.sample == "single") {
        return(row.tokenize('\t')[2])
    } else {
        return(row.tokenize('\t')[4])
    }
}

//  Grab the set of sample IDs from the manifest
sample_ids = file("${params.input}/samples.manifest")
    .readLines()
    .collect{ get_ids(it) }
    .flatten()
    .unique( false )

//  Deduce sample ID from a filename that includes it
def get_prefix(f, trim_read_num=false, check_ids=true) {
    if (trim_read_num) {
        //  Remove these regardless of position in the string
        blackListAny = ~/_[12]_(un|)trimmed_(summary|fastqc_data)|_trimmed|_untrimmed|_unpaired|_paired/
        
        prefix = f.name.toString().replaceAll("_[12]\\.", ".")
    } else {
        //  Remove these regardless of position in the string
        blackListAny = ~/_(un|)trimmed_(summary|fastqc_data)|_trimmed|_untrimmed|_unpaired|_paired/
        
        prefix = f.name.toString()
    }
    
    prefix = prefix
        .replaceAll("\\.Forward", "_Forward")
        .replaceAll("\\.Reverse", "_Reverse")
        .tokenize('.')[0]
        .replaceAll(blackListAny, "")
    
    //  This function returns the sample ID associated with a filename with
    //  typical usage; if the "ID" we found is not one of the sample IDs from
    //  'samples.manifest', something went wrong!
    if (check_ids && ! sample_ids.contains(prefix)) {
        exit 1, "Error: deduced the sample ID '" + prefix + "' from the file '" + f.name.toString() + "', but this ID is not in 'samples.manifest'. This is likely a bug in SPEAQeasy!"
    }
    
    return(prefix)
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

//  Given a "row" of the 'samples.manifest' file as a string, return the FASTQ
//  files
def get_fastq_names(row) {
    if (params.sample == "single") {
        return(file(row.tokenize('\t')[0]))
    } else {
        return(tuple(file(row.tokenize('\t')[0]), file(row.tokenize('\t')[2])))
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Summary of Defined Variables
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This gets the SHA commit ID of the repository where SPEAQeasy is installed.
// This associates the pipeline run with a precise "version" of SPEAQeasy. Note
// that nextflow provides the "workflow.commitId" variable with this intended
// function- during testing this variable appears to be null.
params.commitId = "git --git-dir=${workflow.projectDir}/.git rev-parse HEAD".execute().text.trim()

// Header log info
def summary_main = [:]
summary_main['SPEAQeasy version'] = params.commitId
summary_main['Working dir']      = workflow.workDir
summary_main['Current home']        = "$HOME"
summary_main['Current user']        = "$USER"
summary_main['Current path']        = "$PWD"
summary_main['Annotation build'] = params.anno_build
summary_main['Annotation dir']       = params.annotation
summary_main['Annotation release'] = params.anno_version
summary_main['Compute coverage'] = params.fullCov || params.coverage
summary_main['Custom anno label'] = params.custom_anno
summary_main['ERCC spike-in'] = params.ercc
summary_main['Experiment name'] = params.experiment
summary_main['Full coverage'] = params.fullCov
summary_main['Input dir']              = params.input
summary_main['Keep unpaired'] = params.keep_unpaired
summary_main['Output dir']        = params.output
summary_main['Prefix'] = params.prefix
summary_main['Reference']          = params.reference
summary_main['Sample']            = params.sample
summary_main['Small test']  = params.small_test
summary_main['Strand']            = params.strand
summary_main['Strand mode'] = params.strand_mode
summary_main['Trim mode'] = params.trim_mode
summary_main['Keep discordant'] = params.unalign
summary_main['Use salmon'] = params.use_salmon
summary_main['Use STAR'] = params.use_star

def summary_args = [:]
summary_args['Num reads for strand inference'] = params.num_reads_infer_strand
summary_args['Wiggletools batch size'] = params.wiggletools_batch_size
summary_args['Variants merge batch size'] = params.variants_merge_batch_size
summary_args['bam2wig arguments'] = params.bam2wig_args
summary_args['bcftools arguments'] = params.bcftools_args
summary_args['FastQC arguments'] = params.fastqc_args
summary_args['featureCounts gene arguments'] = params.feat_counts_gene_args
summary_args['featureCounts exon arguments'] = params.feat_counts_exon_args
summary_args['HISAT2 arguments'] = params.hisat2_args
summary_args['kallisto quant single arguments'] = params.kallisto_quant_single_args
summary_args['kallisto quant paired arguments'] = params.kallisto_quant_paired_args
summary_args['kallisto quant ERCC single arguments'] = params.kallisto_quant_ercc_single_args
summary_args['kallisto quant ERCC paired arguments'] = params.kallisto_quant_ercc_paired_args
summary_args['kallisto index arguments'] = params.kallisto_index_args
summary_args['salmon index arguments'] = params.salmon_index_args
summary_args['salmon quant arguments'] = params.salmon_quant_args
summary_args['samtools arguments'] = params.samtools_args
summary_args['STAR arguments'] = params.star_args
summary_args['Adapter trimming single arguments'] = params.trim_adapter_args_single
summary_args['Adapter trimming paired arguments'] = params.trim_adapter_args_paired
summary_args['Quality trimming arguments'] = params.trim_quality_args
summary_args['RegTools arguments'] = params.regtools_args
summary_args['wigToBigWig arguments'] = params.wigToBigWig_args

log.info "================================================================================"
log.info "    SPEAQeasy: an RNA-seq analysis pipeline from LIBD"
log.info "================================================================================"
log.info "---- Main options:"
log.info summary_main.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "---- Software arguments:"
log.info summary_args.collect { k,v -> "${k.padRight(36)}: $v" }.join("\n")
log.info "================================================================================"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN PIPELINE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Pull assembly fasta (from GENCODE for human/ mouse, or ensembl for rat).
 * Build a "main" fasta of only canonical sequences from the pulled "primary"
 * fasta containing all sequences/contigs.
 */

// If any files are not already downloaded/ prepared: download the primary
// assembly fasta and create a subsetted fasta of "main" (reference) sequences only
process PullAssemblyFasta {

    tag "Downloading Assembly FA File: ${baseName}"
    storeDir "${params.annotation}/reference/${params.reference}/assembly/fa"
        
    output:
        path "${out_fasta}", emit: reference_fasta
        path "*.fa" // to store the primary and main fastas, regardless of which is used

    shell:
        //  Name of the primary assembly fasta after being downloaded and unzipped
        baseName = file("${params.fa_link}").getName() - ".gz"
        
        //  Name the pipeline will use for the primary and main assembly fastas, respectively
        primaryName = "assembly_${params.anno_suffix}.fa".replaceAll("main", "primary")
        mainName = "assembly_${params.anno_suffix}.fa".replaceAll("primary", "main")
        
        //  Name of fasta to use for this pipeline execution instance
        out_fasta = "assembly_${params.anno_suffix}.fa"
        
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Pull and unzip primary assembly fasta
        curl -O "!{params.fa_link}"
        gunzip "!{baseName}.gz"
        mv !{baseName} !{primaryName} # rename for consistency with pipeline naming conventions
        
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
        first_bad_line=$(grep -n ">" !{primaryName} | cut -d : -f 1 | paste -s | cut -f $(($num_chrs + 1)))
        
        #  Make a new file out of all the lines up and not including that
        sed -n "1,$(($first_bad_line - 1))p;${first_bad_line}q" !{primaryName} > !{mainName}
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}


/*
 * Download reference .gtf (from GENCODE for human/ mouse, ensembl for rat)
 */
// Uses "storeDir" to download gtf only when it doesn't exist, and output the cached
// file if it does already exist
process PullGtf {

    tag "Downloading GTF File: ${baseName}"
    storeDir "${params.annotation}/RSeQC/${params.reference}/gtf"

    output:
        path "${out_gtf}", emit: reference_gtf

    shell:
        // Names of gtf file when downloaded + unzipped and after renamed, respectively
        baseName = file("${params.gtf_link}").getName() - ".gz"
        out_gtf = "transcripts_${params.anno_suffix}.gtf"
        '''
        #  Pull, unzip, and rename transcript gtf
        curl -O "!{params.gtf_link}"
        gunzip "!{baseName}.gz"
        mv !{baseName} !{out_gtf}
        '''
}

/*
 * Build aligner index (HISAT2 by default, or optionally STAR)
 */

process BuildStarIndex {
    storeDir "${params.annotation}/reference/${params.reference}/assembly/index/star_${params.anno_suffix}"
    
    input:
        path reference_fasta
        path reference_gtf
        
    output:
        path "index_dir/*", emit: star_index
        
    shell:
        '''
        mkdir index_dir
        
        !{params.star} \
            --runMode genomeGenerate \
            --genomeDir ./index_dir \
            --runThreadN !{task.cpus} \
            --genomeFastaFiles !{reference_fasta} \
            --sjdbGTFfile !{reference_gtf}
        '''
}

process BuildHisatIndex {	
    tag "Building HISAT2 Index: ${hisat_prefix}"
    storeDir "${params.annotation}/reference/${params.reference}/assembly/index"

    input:
        path reference_fasta

    output:
        path "${hisat_prefix}.*", emit: hisat_index

    shell:
        hisat_prefix = "hisat2_assembly_${params.anno_suffix}"
        '''
        !{params.hisat2build} -p !{task.cpus} !{reference_fasta} !{hisat_prefix}
        '''
}


// Build "objects" from the assembly fasta and reference gtf. These objects include
// .rda files with junction and transcript feature information, and a file listing
// the size of each chromosome/ sequence. These objects are for internal use by
// the pipeline.
process BuildAnnotationObjects {

    storeDir "${params.annotation}/junction_txdb"

    input:
        path reference_fasta
        path reference_gtf
        path build_ann_script

    output:
        path "junction_annotation_${params.anno_suffix}.rda", emit: junction_annotation
        path "feature_to_Tx_${params.anno_suffix}.rda", emit: feature_to_tx_gencode
        path "chrom_sizes_${params.anno_suffix}", emit: chr_sizes

    shell:
        '''
        !{params.Rscript} !{build_ann_script} -r !{params.reference} -s !{params.anno_suffix}
        '''
}

/*
 * Download transcript FASTA
 */
			
// Uses "storeDir" to download files only when they don't exist, and output the cached
// files if they do already exist
process PullTranscriptFasta {
            
    tag "Downloading TX FA File: ${baseName}"
    storeDir "${params.annotation}/reference/${params.reference}/transcripts/fa"
    
    input:
        path subset_script

    output:
        path baseName, emit: transcript_fa

    shell:
        //  For human and mouse, only the "main" transcripts FASTA is
        //  available. These means the output FASTA won't differ between
        //  "main" and "primary" runs for these species. For rat, only a
        //  "primary" FASTA is available, so the output will differ
        if (params.reference == 'rat') {
            baseName = "transcripts_${params.anno_suffix}.fa"
        } else {
            baseName = file("${params.tx_fa_link}").getName() - ".gz"
        }

        '''
        curl -o !{baseName}.gz !{params.tx_fa_link}
        gunzip !{baseName}.gz

        #   For rat, the only transcripts FASTA available includes "primary"
        #   transcripts. Subset if the user selects "main" build
        if [[ (!{params.reference} == 'rat') && (!{params.anno_build} == "main") ]]; then
            !{params.Rscript} !{subset_script}
        fi
        '''
}

/*
 * Build transcript index for Salmon or Kallisto
 */

// Uses "storeDir" to build the index only if the built file is not present; outputs
// this cached file otherwise
process BuildSalmonIndex {

    tag "Building Salmon Index: salmon_index_${anno_suffix}"
    storeDir "${params.annotation}/reference/${params.reference}/transcripts/salmon/${anno_suffix}"

    input:
        path transcript_fa

    output:
        path "salmon_index_${anno_suffix}/*", emit: salmon_index
        path "build_salmon_index_${anno_suffix}.log"

    shell:
        if (params.reference == "rat") {
            gencode_flag = ""
        } else {
            gencode_flag = "--gencode"
        }

        // Use of main/primary doesn't affect transcripts for human and mouse
        if ((params.custom_anno == "") || (params.reference == "rat")) {
            anno_suffix = params.anno_suffix - "_${params.anno_build}"
        } else {
            anno_suffix = params.anno_suffix
        }
        '''
        !{params.salmon} index \
            -t !{transcript_fa} \
            -i salmon_index_!{anno_suffix} \
            -p !{task.cpus} \
            !{gencode_flag} \
            !{params.salmon_index_args}
        
        cp .command.log build_salmon_index_!{anno_suffix}.log
        '''
}

process BuildKallistoIndex {

    storeDir "${params.annotation}/reference/${params.reference}/transcripts/kallisto"
    
    input:
        path transcript_fa

    output:
        path "kallisto_index_${anno_suffix}", emit: kallisto_index
        path "build_kallisto_index_${anno_suffix}.log"

    shell:
        // Use of main/primary doesn't affect transcripts for human and mouse
        if ((params.custom_anno == "") || (params.reference == "rat")) {
            anno_suffix = params.anno_suffix - "_${params.anno_build}"
        } else {
            anno_suffix = params.anno_suffix
        }

        '''
        !{params.kallisto} index \
            -i kallisto_index_!{anno_suffix} \
            !{params.kallisto_index_args} \
            !{transcript_fa}
        cp .command.log build_kallisto_index_!{anno_suffix}.log
        '''
}


//  Merge FASTQ files as necessary, rename files based on sample ids provided
//  in the manifest, and create a new manifest for internal use based on these
//  changes

process PreprocessInputs {

    publishDir "${params.output}/preprocessing", mode:'copy', pattern:'*.log'

    input:
        path original_manifest
        path merge_script
        path raw_fastqs

    output:
        path "*.f*q*", includeInputs: true, emit: merged_inputs_flat
        path "samples_processed.manifest", emit: strandness_manifest
        path "preprocess_inputs.log"

    shell:
        if (params.sample == "paired") {
            paired_arg = "TRUE"
        } else {
            paired_arg = "FALSE"
        }
        '''
        !{params.Rscript} !{merge_script} -p !{paired_arg}
        cp .command.log preprocess_inputs.log
        '''
}

/*
 * Perform FastQC on the untrimmed reads, as an initial quality gauge
 */

process QualityUntrimmed {

    tag "$prefix"
    publishDir "${params.output}/fastQC/untrimmed", mode:'copy', pattern:'*_fastqc'

    input:
        tuple val(prefix), path(fastqc_untrimmed_input)

    output:
        path "${prefix}*_fastqc"
        tuple val(prefix), path("*_summary.txt"), emit: quality_reports
        path "*_fastqc_data.txt", emit: count_objects_quality_metrics_untrimmed

    shell:
        if (params.sample == "single") {
            copy_command = "cp ${prefix}_fastqc/summary.txt ${prefix}_untrimmed_summary.txt"
            data_command = "cp ${prefix}_fastqc/fastqc_data.txt ${prefix}_untrimmed_fastqc_data.txt"
        } else {
            copy_command = "cp ${prefix}_1_fastqc/summary.txt ${prefix}_1_untrimmed_summary.txt && cp ${prefix}_2_fastqc/summary.txt ${prefix}_2_untrimmed_summary.txt"
            data_command = "cp ${prefix}_1_fastqc/fastqc_data.txt ${prefix}_1_untrimmed_fastqc_data.txt && cp ${prefix}_2_fastqc/fastqc_data.txt ${prefix}_2_untrimmed_fastqc_data.txt"
        }
        '''
        !{params.fastqc} -t !{task.cpus} --extract !{params.fastqc_args} !{fastqc_untrimmed_input}
        !{copy_command}
        !{data_command}
        '''
} 

// Perform pseudoaligment by sample on a subset of reads, and see which strandness
// assumption passed to kallisto results in the largest number of alignments. This
// determines strandness for each sample, used by all remaining pipeline steps which
// require this information.
process InferStrandness {

    tag "${prefix}"
    publishDir "${params.output}/infer_strandness", mode:'copy'
    
    input:
        path infer_strand_R
        path infer_strand_sh
        path kallisto_index
        // The FastQC summary is input here just to ensure FastQC runs before
        // potential strand-inference errors, where FastQC can provide useful
        // insight about potential problems
        tuple val(prefix), path(fastqc_summary), path(fq_file)
        
    output:
        path "${prefix}_infer_strand.log"
        path "${prefix}_strandness_pattern.txt", emit: strandness_patterns
        
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
            !{params.Rscript} \
            !{params.strand_mode} \
            !{prefix}
        
        cp .command.log !{prefix}_infer_strand.log
        '''
}

// Attach strandness information from the InferStrandness process to a copy
// of the user-provided samples.manifest file, for internal use by the
// pipeline.
process CompleteManifest {

    publishDir "${params.output}/infer_strandness", mode:'copy'
    
    input:
        path strandness_files
        path strandness_manifest
        path manifest_script
        
    output:
        path "samples_complete.manifest", emit: complete_manifest
        
    shell:
        '''
        !{params.Rscript} !{manifest_script}
        '''
}

/*
 * Run the ERCC process if the --ercc flag is specified
 */

process ERCC {
	
    tag "$prefix"
    publishDir "${params.output}/ERCC/${prefix}",'mode':'copy'

    input:
        path ercc_index
        tuple val(prefix), path(ercc_input)
        path complete_manifest

    output:
        path "${prefix}_ercc_abundance.tsv", emit: ercc_abundances
        path "ercc_${prefix}.log"

    shell:
        if (params.sample == "single") {
            kallisto_flags = params.kallisto_quant_ercc_single_args
        } else {
            kallisto_flags = params.kallisto_quant_ercc_paired_args
        }
        '''
        ( set -o posix ; set ) > bash_vars.txt

        #  Find this sample's strandness and determine kallisto flags to set
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            kallisto_strand=" --fr-stranded"
        elif [ $strand == 'reverse' ]; then
            kallisto_strand=" --rf-stranded"
        else
            kallisto_strand=""
        fi

        #  Quantify ERCC, and don't throw an error if 0 counts are found
        !{params.kallisto} quant \
            -i !{ercc_index} \
            -t !{task.cpus} \
            !{kallisto_flags} \
            -o . \
            $kallisto_strand \
            !{ercc_input} \
            || true

        cp abundance.tsv !{prefix}_ercc_abundance.tsv
        cp .command.log ercc_!{prefix}.log

        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep '>' | cut -d ' ' -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

/*
 * Trim samples (dependent on user-chosen settings)
 */
 
process Trimming {

    tag "$fq_prefix"
    publishDir "${params.output}/trimming", mode:'copy', pattern:'*_trimmed*.f*q{.gz,}'

    input:
        tuple val(fq_prefix), path(fq_summary), path(fq_file)

    output:
        tuple val(fq_prefix), path("${fq_prefix}_trimmed*.f*q{.gz,}"), optional: true, emit: trimmed_fastqc_inputs_raw
        tuple val(fq_prefix), path("${fq_prefix}*.f*q{.gz,}"), emit: trimming_outputs

    shell:
        file_ext = get_file_ext(fq_file[0])
        if (params.sample == "single") {
            output_option = "${fq_prefix}_trimmed.fastq"
            trim_mode = "SE"
            adapter_fa_temp = params.adapter_fasta_single
            trim_clip = params.trim_adapter_args_single
        } else {
            output_option = "${fq_prefix}_trimmed_paired_1.fastq ${fq_prefix}_unpaired_1.fastq ${fq_prefix}_trimmed_paired_2.fastq ${fq_prefix}_unpaired_2.fastq"
            trim_mode = "PE"
            adapter_fa_temp = params.adapter_fasta_paired
            trim_clip = params.trim_adapter_args_paired
        }
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Determine whether to trim the FASTQ file(s). This is ultimately
        #  controlled by the '--trim_mode' command flag.
        if [ "!{params.trim_mode}" == "force" ]; then
            do_trim=true
        elif [ "!{params.trim_mode}" == "skip" ]; then
            do_trim=false
        elif [ "!{params.sample}" == "single" ]; then
            #  Then '--trim_mode "adaptive"' was selected, and data is single-end
            if [ $(grep "Adapter Content" !{fq_summary} | cut -f 1)  == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        else
            #  Then '--trim_mode "adaptive"' was selected, and data is paired-end
            result1=$(grep "Adapter Content" !{fq_prefix}_1_untrimmed_summary.txt | cut -c1-4)
            result2=$(grep "Adapter Content" !{fq_prefix}_2_untrimmed_summary.txt | cut -c1-4)
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
            
            #  Now determine the arguments to pass to trimmomatic regarding
            #  adapter trimming. The flexibility here allows the user to 
            #  potentially not perform adapter trimming.
            if [ "!{trim_clip}" == "" ]; then
                adapter_trim_settings=""
            else
                adapter_trim_settings="ILLUMINACLIP:$adapter_fa:!{trim_clip}"
            fi

            java -Xmx512M \
                -jar $trim_jar \
                !{trim_mode} \
                -threads !{task.cpus} \
                -phred33 \
                *.f*q* \
                !{output_option} \
                $adapter_trim_settings \
                !{params.trim_quality_args}
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
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

/*
 * Run FastQC Quality Check on Trimmed Files
 */

process QualityTrimmed {

    tag "$prefix"
    publishDir "${params.output}/fastQC/trimmed",'mode':'copy', pattern:'*_fastqc'

    input:
        tuple val(prefix), path(fastqc_trimmed_input)

    output:
        path "${prefix}*_fastqc"
        path "${prefix}*_trimmed_*.txt", emit: trimmed_fastqc_outs

    shell:
        if (params.sample == "single") {
            copy_command = "cp ${prefix}_trimmed_fastqc/summary.txt ${prefix}_trimmed_summary.txt"
            data_command = "cp ${prefix}_trimmed_fastqc/fastqc_data.txt ${prefix}_trimmed_fastqc_data.txt"
        } else {
            copy_command = "cp ${prefix}_trimmed_paired_1_fastqc/summary.txt ${prefix}_1_trimmed_summary.txt && cp ${prefix}_trimmed_paired_2_fastqc/summary.txt ${prefix}_2_trimmed_summary.txt"
            data_command = "cp ${prefix}_trimmed_paired_1_fastqc/fastqc_data.txt ${prefix}_1_trimmed_fastqc_data.txt && cp ${prefix}_trimmed_paired_2_fastqc/fastqc_data.txt ${prefix}_2_trimmed_fastqc_data.txt"
        }
        '''
        !{params.fastqc} -t !{task.cpus} --extract !{params.fastqc_args} !{fastqc_trimmed_input}
        !{copy_command}
        !{data_command}
        '''
}

/*
 * Alignment (HISAT2 by default, or otherwise STAR)
 */

process AlignStar {
    tag "$prefix"
    publishDir "${params.output}/alignment", mode:'copy'
    
    input:
        path star_index
        tuple val(prefix), path(single_star_input)
        
    output:
        tuple val(prefix), path("*.bam"), emit: alignment_output
        path "${prefix}_unmapped_*.fastq", optional: true
        path "*_STAR_alignment.log", emit: alignment_summaries
        
    shell:            
        if (params.sample == "paired" && params.unalign) {
            unaligned_args = "--outReadsUnmapped Fastx"
        } else {
            unaligned_args = ""
        }
        '''
        #  Recreate the index directory (nextflow cannot handle directories
        #  as channel input). This step is likely sensitive to the STAR
        #  version used
        mkdir index_dir
        mv *.txt index_dir/
        mv *.tab index_dir/
        mv SA index_dir/
        mv SAindex index_dir/
        mv Genome index_dir/
        
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Determine file names for input FASTQs
        if [ !{params.sample} == "paired" ]; then
            fastq_files="!{prefix}*trimmed*_1.f*q* !{prefix}*trimmed*_2.f*q*"
        else
            fastq_files="*.f*q*"
        fi
        
        #  Tell STAR to decompress files if needed 
        if [ "$(echo $fastq_files | rev | cut -d "." -f 1 | rev)" == "gz" ]; then
            compress_args='--readFilesCommand gunzip -c'
        else
            compress_args=''
        fi
        
        #  Perform alignment
        !{params.star} \
            --genomeDir ./index_dir \
            --runThreadN !{task.cpus} \
            --readFilesIn ${fastq_files} \
            --outSAMtype BAM Unsorted \
            ${compress_args} \
            !{unaligned_args} \
            !{params.star_args}
            
        #  Adjust file names (make unique)
        mv Log.final.out !{prefix}_STAR_alignment.log
        mv Aligned.out.bam !{prefix}.bam
        if [ -f Unmapped.out.mate1 ]; then
            mv Unmapped.out.mate1 !{prefix}_unmapped_mate1.fastq
            mv Unmapped.out.mate2 !{prefix}_unmapped_mate2.fastq
        fi
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}
        
process SingleEndHisat {

    tag "$prefix"
    publishDir "${params.output}/alignment",mode:'copy'

    input:
        path hisat_index
        path complete_manifest
        tuple val(prefix), path(single_hisat_input)

    output:
        tuple val(prefix), path("${prefix}.bam"), emit: alignment_output
        path "*_align_summary.txt", emit: alignment_summaries

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
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
            ${hisat_strand} \
            !{params.hisat2_args} \
            2> !{prefix}_align_summary.txt \
        | !{params.samtools} view -b -F 4 -o !{prefix}.bam
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

process PairedEndHisat {

    tag "$prefix"
    publishDir "${params.output}/alignment",'mode':'copy'

    input:
        path hisat_index
        path complete_manifest
        tuple val(prefix), path(fq_files)

    output:
        tuple val(prefix), path("${prefix}.bam"), emit: alignment_output
        path "*_unpaired*.fastq" optional true
        path "*_align_summary.txt", emit: alignment_summaries

    shell:
        if (params.unalign) {
            unaligned_opt = "--un-conc ${prefix}_discordant.fastq"
        } else {
            unaligned_opt = ""
        }
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ ${strand} == "unstranded" ]; then
            hisat_strand=""
        elif [ ${strand} == "forward" ]; then
            hisat_strand="--rna-strandness FR"
        else
            hisat_strand="--rna-strandness RF"
        fi
        
        #  If this sample had unpaired trimming outputs, include them
        if [ "!{params.keep_unpaired}" == "true" ]; then
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
            ${unpaired_opt} \
            ${hisat_strand} \
            !{params.hisat2_args} \
            !{unaligned_opt} \
            2> !{prefix}_align_summary.txt \
        | !{params.samtools} view -b -F 4 -o !{prefix}.bam
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

/*
 * Sort and index the aligned BAM
 */

process BamSort {
	
    tag "$prefix"
    publishDir "${params.output}/alignment/bam_sort",'mode':'copy'

    input:
        tuple val(prefix), path(input_bam)

    output:
        tuple val(prefix), path("${prefix}_sorted.bam"), path("${prefix}*_sorted.bam.bai"), emit: sorted_bams

    shell:
        '''
        !{params.samtools} sort -T temporary -l 6 --no-PG -@ !{task.cpus} !{input_bam} -o !{prefix}_sorted.bam
        !{params.samtools} index !{prefix}_sorted.bam
        '''
}

/*
 * Quantify genes and exons with FeatureCounts
 */

process FeatureCounts {

    tag "$prefix"
    publishDir "${params.output}/counts",'mode':'copy'

    input:
        tuple val(prefix), path(bam), path(bam_index)
        path complete_manifest
        path reference_gtf

    output:
        path "*.{log,summary}"
        path "*.counts*", emit: sample_counts

    script:
        if (params.sample == "single") {
            sample_option = ""
        } else {
            sample_option = "-p"
        }
        feature_out = "${prefix}_${params.anno_suffix}"
 
        """
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=\$(cat samples_complete.manifest | grep " ${prefix} " | awk -F ' ' '{print \$NF}')
        if [ \$strand == 'forward' ]; then
            feature_strand=1
        elif [ \$strand == 'reverse' ]; then
            feature_strand=2
        else
            feature_strand=0
        fi
        
        # Genes
        ${params.featureCounts} \
            -s \$feature_strand \
            $sample_option \
            ${params.feat_counts_gene_args} \
            -T $task.cpus \
            -a $reference_gtf \
            -o ${feature_out}_Genes.counts \
            $bam
        
        # Exons
        ${params.featureCounts} \
            -s \$feature_strand \
            $sample_option \
            ${params.feat_counts_exon_args} \
            -T $task.cpus \
            -a $reference_gtf \
            -o ${feature_out}_Exons.counts \
            $bam

        cp .command.log feature_counts_${prefix}.log
        
        temp=\$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "\$temp" > bash_vars.txt
        """
}

/*
 * Primary Alignments
 */

process PrimaryAlignments {

    tag "$alignment_prefix"
    publishDir "${params.output}/counts/junction/primary_alignments",'mode':'copy'

    input:
        tuple val(alignment_prefix), path(alignment_bam), path(alignment_index)

    output:
        tuple val("${alignment_prefix}"), path("${alignment_prefix}.bam"), path("${alignment_prefix}.bam.bai"), emit: primary_alignments

    script:
        """
        ${params.samtools} view -@ $task.cpus -bh -F 0x100 $alignment_bam > ${alignment_prefix}.bam
        ${params.samtools} index ${alignment_prefix}.bam
        """
}

/*
 * Quantify Junctions
 */

process Junctions {

    tag "$prefix"
    publishDir "${params.output}/counts/junction",'mode':'copy'

    input:
        tuple val(prefix), path(alignment_bam), path(alignment_index)
        path complete_manifest

    output:
        path "*.count", emit: regtools_outputs

    shell:
        outcount = "${prefix}_junctions_primaryOnly_regtools.count"
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            strand_integer=2
        elif [ $strand == 'reverse' ]; then
            strand_integer=1
        else
            strand_integer=0
        fi
        
        !{params.regtools} junctions extract !{params.regtools_args} -s ${strand_integer} -c !{outcount} !{alignment_bam}
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}


// /*
//  * Transcript quantification
//  */

process TxQuantSalmon {

    tag "$prefix"
    publishDir "${params.output}/salmon_tx/${prefix}",mode:'copy'

    input:
        path salmon_index
        path complete_manifest
        tuple val(prefix), path(salmon_inputs)

    output:
        path "${prefix}/*"
        path "${prefix}_quant.sf", emit: tx_quants
        path "tx_quant_${prefix}.log"

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
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
            -i $PWD \
            -p !{task.cpus} \
            -l ${strand_flag} \
            ${sample_flag} \
            -o !{prefix} \
            !{params.salmon_quant_args}

        cp !{prefix}/quant.sf !{prefix}_quant.sf
        cp .command.log tx_quant_!{prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

process TxQuantKallisto {

    tag "$prefix"
    publishDir "${params.output}/kallisto_tx/${prefix}", mode:'copy'

    input:
        path kallisto_index
        path complete_manifest
        tuple val(prefix), path(fastqs)

    output:
        path "abundance.h5"
        path "run_info.json"
        path "tx_quant_${prefix}.log"
        path "${prefix}_abundance.tsv", emit: tx_quants
        
    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
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
        
            !{params.kallisto} quant \
                -t !{task.cpus} \
                $strand_flag \
                -i !{kallisto_index} \
                -o . \
                !{params.kallisto_quant_paired_args} \
                $fq1 $fq2
        else
            fq=$(ls *.f*q*)
            
            !{params.kallisto} quant \
                -t !{task.cpus} \
                !{params.kallisto_quant_single_args} \
                $strand_flag \
                -i !{kallisto_index} \
                -o . \
                $fq
        fi
        
        mv abundance.tsv !{prefix}_abundance.tsv
        cp .command.log tx_quant_!{prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

/*
 * Create RangedSummarizedExperiment outputs
*/

process CountObjects {

    publishDir "${params.output}/count_objects",'mode':'copy'

    input:
        // Reference/ annotation-related files
        path ercc_reference
        path exon_maps_by_coord_hg38
        path reference_gtf
        path feature_to_tx_gencode
        path junction_annotation
        // FastQC summaries before and after trimming
        path count_objects_quality_metrics_untrimmed
        path quality_reports
        path trimmed_fastqc_outs
        // Summaries and BAMs from alignment
        path alignment_summaries
        path bams_and_indices
        // Gene + exon counts
        path sample_counts
        // Junction counts
        path regtools_outputs
        // Transcript counts
        path tx_quants
        // File containing list of transcripts, when using --qsva
        path qsva_tx_list
        // ERCC spike-in counts (if applicable)
        path ercc_abundances
        // Manifest with strandness info
        path complete_manifest
        // R script for creating the RSE objects
        path counts_script


    output:
        path "*.pdf", optional: true
        path "read_and_alignment_metrics_*.csv"
        path "*.Rdata"
        path "*.rda"
        path "counts.log"

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
        # Write 'params' to CSV, where it can be read in (in R) and used to
        # record SPEAQeasy settings in each RSE's metadata
        echo "!{params}" | sed 's|, |\\n|g' | tr -d '[]' | sed 's|:|,|' > params.csv
        
        if [[ "!{params.qsva}" == "" ]]; then
            qsva_arg=""
        else
            qsva_arg="-q $(basename !{params.qsva})"
        fi
        
        !{params.Rscript} !{counts_script} \
            -o !{params.reference} \
            -e !{params.experiment} \
            -p "!{params.prefix}" \
            -l !{counts_pe} \
            -c !{params.ercc} \
            -t !{task.cpus} \
            !{counts_strand} \
            -n !{params.use_salmon} \
            -r !{params.use_star} \
            -u !{params.output} \
            -a !{params.anno_suffix} \
            ${qsva_arg}

        cp .command.log counts.log
        '''
}

process VariantCalls {
    tag "$prefix"
    publishDir "${params.output}/variant_calls",'mode':'copy'
    
    input:
        tuple val(prefix), path(sorted_bams_file), path(variant_calls_bai)
        path snv_bed
        path reference_fasta
    
    output:
        path "${prefix}.vcf.gz", emit: compressed_variant_calls
        path "${prefix}.vcf.gz.tbi", emit: compressed_variant_calls_tbi
    
    shell:
        '''
        !{params.samtools} mpileup \
            -l !{snv_bed} \
            !{params.samtools_args} \
            -u \
            -f !{reference_fasta} \
            !{sorted_bams_file} \
            -o !{prefix}_tmp.vcf
        
        !{params.bcftools} call \
            !{params.bcftools_args} \
            !{prefix}_tmp.vcf \
            > !{prefix}.vcf.gz
        
        !{params.tabix} -p vcf !{prefix}.vcf.gz
        '''
}


/*
 * Merge Variant Calls
 */

process VariantsMerge {

    tag "Multi-sample vcf creation"
    publishDir "${params.output}/merged_variants",'mode':'copy'

    input:
        path chunk_apply_script
        path collected_variants
        path collected_variants_tbi

    output:
        path "mergedVariants.vcf.gz"
        path "variants_merge.log"

    shell:
        '''
        file_regex='.*\\.vcf\\.gz$'
        command="!{params.bcftools} merge [files] | !{params.bgzip} -c > temp[index].vcf.gz; !{params.tabix} -p vcf temp[index].vcf.gz"

        #   Break the 'bcftools merge' command into chunks (necessary for large
        #   datasets where the number of open file handles may be exceeded)
        bash !{chunk_apply_script} !{params.variants_merge_batch_size} "$file_regex" "$command"
        
        file_list=$(ls -1 | grep -E "$file_regex")
        num_files=$(echo "$file_list" | wc -l)
        if [[ $num_files -gt 1 ]]; then
            #   If more than one file is left, merge everything to produce the
            #   final output
            echo "Performing a final merge of all VCF files..."
            !{params.bcftools} merge $file_list | !{params.bgzip} -c > mergedVariants.vcf.gz
        else
            #   There's only one file, so rename it
            mv $file_list mergedVariants.vcf.gz
        fi

        cp .command.log variants_merge.log
        '''
}

process Coverage {

    tag "$prefix"
    publishDir "${params.output}/coverage/wigs", mode:'copy'

    input:
        path complete_manifest
        tuple val(prefix), path(sorted_bam), path(sorted_bam_index)
        path chr_sizes

    output:
        tuple val(prefix), path("${prefix}*.wig"), emit: wig_files
        path "bam2wig_${prefix}.log"

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            if [ !{params.sample} == "paired" ]; then
                strand_flag='-d 1++,1--,2+-,2-+'
            else
                strand_flag='-d ++,--'
            fi
        elif [ $strand == 'reverse' ]; then
            if [ !{params.sample} == "paired" ]; then
                strand_flag='-d 1+-,1-+,2++,2--'
            else
                strand_flag='-d +-,-+'
            fi
        else
            strand_flag=""
        fi

        python3 $(which bam2wig.py) \
            -s !{chr_sizes} \
            -i !{sorted_bam} \
            !{params.bam2wig_args} \
            -o !{prefix} \
            $strand_flag

        cp .command.log bam2wig_!{prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2-|| true)
        echo "$temp" > bash_vars.txt
        '''
}

// Convert wig files to bigWig format    
process WigToBigWig {

    tag "$prefix"
    publishDir "${params.output}/coverage/bigWigs", mode: 'copy'
    
    input:
        tuple val(prefix), path(wig_file)
        path chr_sizes
    
    output:
        path "*.bw", emit: coverage_bigwigs
    
    shell:
        '''
        !{params.wigToBigWig} \
            !{params.wigToBigWig_args} \
            !{wig_file} \
            !{chr_sizes} \
            !{prefix}.bw
        '''
}

// Compute mean coverage across samples by strand
process MeanCoverage {

    tag "Strand: ${read_type}"
    publishDir "${params.output}/coverage/mean", mode: 'copy'

    input:
        tuple val(read_type), path(mean_coverage_bigwig)
        path chr_sizes
        path chunk_apply_script

    output:
        tuple val(read_type), path("mean*.bw"), emit: mean_bigwigs
        path "mean_coverage_${read_type}.log"

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        if [ !{read_type} == "unstranded" ] ; then
            outwig="mean"
        elif [ !{read_type} == "forward" ]; then
            outwig="mean.forward"
        else
            outwig="mean.reverse"
        fi
        
        precision=8
        
        file_regex='.*\\.(bw|wig)$'
        command="!{params.wiggletools} write temp_wig[index].wig sum [files]"

        #   To compute a mean, we'll sum all files then multiply by the
        #   scale factor (reciprocal of number of files)
        num_files=$(ls -1 | grep -E "$file_regex" | wc -l)
        scale_factor=$(!{params.bc} <<< "scale=${precision};1/$num_files")
        echo "Scale factor is $scale_factor."
        
        #   Break the summation command into chunks (necessary for large
        #   datasets where the number of open file handles may be exceeded)
        bash !{chunk_apply_script} !{params.wiggletools_batch_size} "$file_regex" "$command"

        #   Note that even if one file is left, this command is necessary
        #   because of the scaling (mean instead of sum)
        echo "Computing a mean from the final files..."
        file_list=$(ls -1 | grep -E "$file_regex")
        !{params.wiggletools} write $outwig.wig scale $scale_factor sum $file_list
        echo "Done."
        
        !{params.wigToBigWig} $outwig.wig !{chr_sizes} $outwig.bw
        
        cp .command.log mean_coverage_!{read_type}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}


process CoverageObjects {
    
    tag "Strand: ${read_type}"
    publishDir "${params.output}/coverage_objects", mode: 'copy'

    input:
        path fullCov_script
        path complete_manifest
        tuple val(read_type), path(bigwigs)

    output:
        file "*"

    shell:
        '''
        !{params.Rscript} !{fullCov_script} \
            -o !{params.reference} \
            -e !{params.experiment} \
            -c !{task.cpus} \
            -s !{read_type}
        
        cp .command.log coverage_objects.log
        '''
}

process ExpressedRegions {
    
    tag "Strand: ${read_type}"
    publishDir "${params.output}/expressed_regions", mode:'copy'
    
    input:
        path expressed_regions_script
        path chr_sizes
        tuple val(read_type), path(mean_bigwigs)
    
    output:
        file "*"
    
    shell:
        '''
        for meanfile in ./mean*.bw; do
            !{params.Rscript} !{expressed_regions_script} \
                -m ${meanfile} \
                -o . \
                -i !{chr_sizes} \
                -c !{task.cpus}
        done
        
        cp .command.log expressed_regions_!{read_type}.log
        '''
}


workflow {
    if (params.reference == "hg38" && params.anno_version == "25") {
        // Our extra exon annotation if user defaults to gencode release 25
        exon_maps_by_coord_hg38 = Channel.fromPath("${params.annotation}/junction_txdb/exonMaps_by_coord_hg38_gencode_v25.rda")
    } else {
        // Nextflow requires a valid path of length 1: here we just use a random script
        exon_maps_by_coord_hg38 = Channel.fromPath("${workflow.projectDir}/scripts/style_code.R")
    }

    // When using custom annotation, grab the user-provided reference FASTA.
    // Otherwise, run PullAssemblyFasta to pull it from online
    if (params.custom_anno != "") {
        reference_fasta = Channel.fromPath("${params.annotation}/*assembly*.fa*")
            .ifEmpty{ error "Cannot find assembly fasta in annotation directory (and --custom_anno was specified)" }
            .first()  // This proves to nextflow that the channel will always hold one value/file
            .collect()
        reference_gtf = Channel.fromPath("${params.annotation}/*.gtf")
            .ifEmpty{ error "Cannot find reference gtf in annotation directory (and --custom_anno was specified)" }
            .first()
            .collect()
        transcript_fa = Channel.fromPath("${params.annotation}/*transcripts*.fa*")
            .ifEmpty{ error "Cannot find transcripts fasta in annotation directory (and --custom_anno was specified)" }
            .first()
            .collect()
        ercc_index = Channel.fromPath("${params.annotation}/*.idx").collect()

        // Variant calling is only enabled for human
        if (params.reference_type == "human") {
            snv_bed = Channel.fromPath("${params.annotation}/*.bed").collect()
        }
    } else {
        PullAssemblyFasta()
        PullGtf()
        PullTranscriptFasta(
            Channel.fromPath("${workflow.projectDir}/scripts/subset_rat_fasta.R")
        )
        reference_fasta = PullAssemblyFasta.out.reference_fasta.collect()
        reference_gtf = PullGtf.out.reference_gtf.collect()
        transcript_fa = PullTranscriptFasta.out.transcript_fa.collect()
        ercc_index = Channel.fromPath("${params.annotation}/ERCC/ERCC92.idx").collect()

        // Note that variant calling is only performed for human
        snv_bed = Channel.fromPath("${params.annotation}/Genotyping/common_missense_SNVs_${params.reference}.bed").collect()
    }

    BuildAnnotationObjects(
        reference_fasta,
        reference_gtf,
        Channel.fromPath("${workflow.projectDir}/scripts/build_annotation_objects.R")
    )

    if (params.use_star) {
        BuildStarIndex(reference_fasta, reference_gtf)
    } else { // HISAT2 is used as the aligner
        BuildHisatIndex(reference_fasta)
    }

    if (params.use_salmon) {
        BuildSalmonIndex(transcript_fa)
    }
    BuildKallistoIndex(transcript_fa) // Needed either way, because of InferStrandness

    // Extract FASTQ file paths from the manifest and place in a channel to pass to
    // PreprocessInputs
    original_manifest = Channel.fromPath("${params.input}/samples.manifest")
    original_manifest
        .splitText()
        .map{ row -> get_fastq_names(row) }
        .flatten()
        .collect()
        .set{ raw_fastqs }
    PreprocessInputs(
        original_manifest,
        Channel.fromPath("${workflow.projectDir}/scripts/preprocess_inputs.R"),
        raw_fastqs
    )

    PreprocessInputs.out.merged_inputs_flat
        .flatten()
        .map{file -> tuple(get_prefix(file, params.sample == "paired"), file) }
        .groupTuple()
        .ifEmpty{ error "Input fastq files (after any merging) are missing from the channel"}
        .set{ untrimmed_fastq_files }

    QualityUntrimmed(untrimmed_fastq_files)
    QualityUntrimmed.out.quality_reports
        .join(untrimmed_fastq_files)
        .set{ trimming_inputs }
    
    InferStrandness(
        // The ".collect" is a workaround to make a channel of paths become a
        // value channel, and thus allow re-use across multiple samples
        Channel.fromPath("${workflow.projectDir}/scripts/infer_strand.R").collect(),
        Channel.fromPath("${workflow.projectDir}/scripts/infer_strand.sh").collect(),
        BuildKallistoIndex.out.kallisto_index.collect(),
        trimming_inputs
    )

    CompleteManifest(
        InferStrandness.out.strandness_patterns.collect(),
        PreprocessInputs.out.strandness_manifest,
        Channel.fromPath("${workflow.projectDir}/scripts/complete_manifest.R")
    )

    if (params.ercc) {
        ERCC(
            ercc_index,
            untrimmed_fastq_files,
            CompleteManifest.out.complete_manifest.collect()
        )
        ercc_abundances = ERCC.out.ercc_abundances
    } else {
        // Nextflow requires a valid path of length 1: here we just use a random script
        ercc_abundances = Channel.fromPath("${workflow.projectDir}/scripts/track_runs.sh")
    }

    Trimming(trimming_inputs)
    QualityTrimmed(Trimming.out.trimmed_fastqc_inputs_raw)

    if (params.use_star) {
        AlignStar(
            BuildStarIndex.out.star_index.collect(),
            Trimming.out.trimming_outputs
        )
        alignment_output_temp = AlignStar.out.alignment_output
        alignment_summaries = AlignStar.out.alignment_summaries
    } else {
        // Then aligning with HISAT2
        if (params.sample == "single") {
            SingleEndHisat(
                BuildHisatIndex.out.hisat_index.collect(),
                CompleteManifest.out.complete_manifest.collect(),
                Trimming.out.trimming_outputs
            )
            alignment_output_temp = SingleEndHisat.out.alignment_output
            alignment_summaries = SingleEndHisat.out.alignment_summaries
        } else {
            PairedEndHisat(
                BuildHisatIndex.out.hisat_index.collect(),
                CompleteManifest.out.complete_manifest.collect(),
                Trimming.out.trimming_outputs
            )
            alignment_output_temp = PairedEndHisat.out.alignment_output
            alignment_summaries = PairedEndHisat.out.alignment_summaries
        }
    }

    BamSort(alignment_output_temp)
    FeatureCounts(
        BamSort.out.sorted_bams,
        CompleteManifest.out.complete_manifest.collect(),
        reference_gtf
    )

    PrimaryAlignments(BamSort.out.sorted_bams)
    Junctions(
        PrimaryAlignments.out.primary_alignments,
        CompleteManifest.out.complete_manifest.collect()
    )

    if (params.use_salmon) {
        TxQuantSalmon(
            BuildSalmonIndex.out.salmon_index.collect(),
            CompleteManifest.out.complete_manifest.collect(),
            untrimmed_fastq_files
        )
        tx_quants = TxQuantSalmon.out.tx_quants
    } else {
        TxQuantKallisto(
            BuildKallistoIndex.out.kallisto_index.collect(),
            CompleteManifest.out.complete_manifest.collect(),
            untrimmed_fastq_files
        )
        tx_quants = TxQuantKallisto.out.tx_quants
    }

    if (params.qsva == "") {
        // Nextflow requires a valid path of length 1: here we just use a random script
        qsva_tx_list = Channel.fromPath("${workflow.projectDir}/scripts/chunk_apply.sh")
    } else {
        qsva_tx_list = Channel.fromPath(params.qsva)
    }

    CountObjects(
        // Reference/ annotation-related files
        Channel.fromPath("${params.annotation}/ERCC/ercc_actual_conc.txt"),
        exon_maps_by_coord_hg38,
        reference_gtf,
        BuildAnnotationObjects.out.feature_to_tx_gencode,
        BuildAnnotationObjects.out.junction_annotation,
        // FastQC summaries before and after trimming
        QualityUntrimmed.out.count_objects_quality_metrics_untrimmed.collect(),
        QualityUntrimmed.out.quality_reports
            // Take just the summaries, not sample IDs
            .map{this_tuple -> this_tuple[1] }
            .collect(),
        QualityTrimmed.out.trimmed_fastqc_outs.collect(),
        // Summaries and BAMs from alignment
        alignment_summaries.collect(),
        BamSort.out.sorted_bams
            // Take just the .bam and .bam.bai files (not the sample ID)
            .map{this_tuple -> tuple(this_tuple[1], this_tuple[2]) }
            .collect(),
        // Gene + exon counts
        FeatureCounts.out.sample_counts.collect(),
        // Junction counts
        Junctions.out.regtools_outputs.collect(),
        // Transcript counts
        tx_quants.collect(),
        // File containing list of transcripts, when using --qsva
        qsva_tx_list,
        // ERCC spike-in counts (if applicable)
        ercc_abundances.collect(),
        // Manifest with strandness info
        CompleteManifest.out.complete_manifest,
        // R script for creating the RSE objects
        Channel.fromPath("${workflow.projectDir}/scripts/create_count_objects.R")   
    )
    if (params.reference_type == "human") {
        VariantCalls(
            BamSort.out.sorted_bams,
            snv_bed,
            reference_fasta
        )
        VariantsMerge(
            Channel.fromPath("${workflow.projectDir}/scripts/chunk_apply.sh"),
            VariantCalls.out.compressed_variant_calls.collect(),
            VariantCalls.out.compressed_variant_calls_tbi.collect()
        )
    }

    //  Compute coverage by sample and also averaged across samples for each
    //  strand
    if (params.coverage || params.fullCov) {
        //  Produce wig files, and convert to bigwig, starting from the
        //  alignment BAMs
        Coverage(
            CompleteManifest.out.complete_manifest,
            BamSort.out.sorted_bams,
            BuildAnnotationObjects.out.chr_sizes.collect()
        )
        WigToBigWig(
            Coverage.out.wig_files,
            BuildAnnotationObjects.out.chr_sizes.collect()
        )

        //  Group by strand and compute a mean across samples
        WigToBigWig.out.coverage_bigwigs
            .map { f -> tuple(get_read_type(f), f) }
            .groupTuple()
            .set{ coverage_bigwigs }
        MeanCoverage(
            coverage_bigwigs,
            BuildAnnotationObjects.out.chr_sizes.collect(),
            Channel.fromPath("${workflow.projectDir}/scripts/chunk_apply.sh").collect()
        )

        if (params.fullCov) {
            CoverageObjects(
                Channel.fromPath("${workflow.projectDir}/scripts/create_fullCov_object.R").collect(),
                CompleteManifest.out.complete_manifest.collect(),
                coverage_bigwigs
            )
        }

        ExpressedRegions(
            Channel.fromPath("${workflow.projectDir}/scripts/find_expressed_regions.R").collect(),
            BuildAnnotationObjects.out.chr_sizes.collect(),
            MeanCoverage.out.mean_bigwigs
        )
    }
}
