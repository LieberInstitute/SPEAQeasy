#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
===================================================================================================================================
                       LIEBER INSTITUTE JAFFE-LAB     R N A - S E Q     A N A L Y S I S     P I P E L I N E  
===================================================================================================================================
 RNA-Seq Multi-Input Analysis Pipeline. Nextflow Version: Started December 2017.
 #### Homepage / Documentation
 https://github.com/LieberInstitute/RNAsp
 #### Authors
 ##### Original Pipeline
 Emily Burke <address@email.com>
 Leonardo Collado-Tores <fellgernon@gmail.com>
 Andrew Jaffee <address@email.com>
 BaDoi Phan <address@email.com>
 ##### Nextflow Version
 Jacob Leonard <leonard.jacob09@gmail.com>
 Israel Aguilar <iaguilaror@gmail.com>
 Violeta Larios <siedracko@gmail.com>
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
    -  C1: Individual Sample Manifest
    -  C2: Sample Manifest
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

    nextflow run main.nf --sample "single" --strand "unstranded" --reference "hg19" --merge --ercc --fullCov -profile standard

    NOTES: The pipeline accepts single or paired end reads. These reads can be stranded or unstranded, and must follow the following naming structure:

    merge: "*_read{1,2}.fastq.gz"
    paired: "*_{1,2}.fastq.gz"

    NOTE: File names can not have "." in the name due to the prefix functions in between process

    nextflow run main.nf {CORE} {OPTIONS}
         {CORE}:
            --sample "single/paired"
              single <- reads files individually "*"
              paired <- reads files paired "*_{1,2}"
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
            --merge <- assumes fastq.gz files require merging "*_read{1,2}.fastq.gz"
	##TODO(iaguilar): Check that THIS IS what --merge means (Docummentation ######)
            --ercc  <- performs ercc quantification
            --fullCov <- performs fullCov R analysis
	    --help <- shows this message
             etc...
	##TODO(iaguilar): Finish the brief help descriptors (Docummentation ######)
    
    Mandatory Options:
    -----------------------------------------------------------------------------------------------------------------------------------
    --sample                      Runs the pipeline for "single" or "paired" end reads
    -----------------------------------------------------------------------------------------------------------------------------------
    --strand                      Runs pipeline for "unstranded, forward, reverse" read types
    -----------------------------------------------------------------------------------------------------------------------------------
    --reference                   Select the desired reference for the run (hg38, hg19, mm10, rn6)
    -----------------------------------------------------------------------------------------------------------------------------------

    Optional Parameters:
    -----------------------------------------------------------------------------------------------------------------------------------
    --experiment                  Name of the experiment being run (ex: "alzheimer"). Defaults to FALSE
    -----------------------------------------------------------------------------------------------------------------------------------
    --prefix                      Defines the prefix of the input files (not used to detect files)
    -----------------------------------------------------------------------------------------------------------------------------------
    --input                       Defines the input folder for the files. Defaults to "./input"
    -----------------------------------------------------------------------------------------------------------------------------------
    --output                      Defines the output folder for the files. Defaults to "./results"
    -----------------------------------------------------------------------------------------------------------------------------------
    --merge                       Flag option if files need to be merged. Defaults to FALSE
    -----------------------------------------------------------------------------------------------------------------------------------
    --unalign                     Give the option to not algin the reads against a reference in HISAT step. Defaults to FALSE 
    -----------------------------------------------------------------------------------------------------------------------------------
    --annotation                  Path to the folder containing pipeline annotations. Defaults to "./Annotations"
    -----------------------------------------------------------------------------------------------------------------------------------
    --indexing                    Path to the base directory containing pipeline indexes. Defaults to --annotation path
    -----------------------------------------------------------------------------------------------------------------------------------
    --genotype                    Path to the folder containing pipeline genotypes. Defaults to "./Genotyping"
	##TODO(iaguilar): Expand genotype description. what are they used for? (Docummentation ######)
    -----------------------------------------------------------------------------------------------------------------------------------
    --email                       Parameter to get a summary e-mail with details of the run sent to you when the workflow exits
	##TODO(iaguilar): Check if this works for failed executions (Docummentation ######)
    -----------------------------------------------------------------------------------------------------------------------------------
    --name                        Name for the pipeline run. If not specified, name is set to experiment
	##TODO(iaguilar): What is this used for? (Docummentation ######)
    -----------------------------------------------------------------------------------------------------------------------------------
    --ercc                        Flag to enable ERCC quantification with Kallisto
    -----------------------------------------------------------------------------------------------------------------------------------
    --k_lm                        Kallisto ERCC Length Mean Value for Single End Reads (defaults to 200)
    -----------------------------------------------------------------------------------------------------------------------------------
    --k_sd                        Kallisto ERCC Standard Deviation Value for Single End Reads (defaults to 30)
    -----------------------------------------------------------------------------------------------------------------------------------
    --fullCov                     Flag to perform full coverage in step 7b
    -----------------------------------------------------------------------------------------------------------------------------------
#################
	##TODO(iaguilar): change this commands to perform quick tests on every posible version of the input data (maybe a quick human test, and a whole species tests)(Docummentation ######)
    --small_test                  Runs the pipeline as a small test run on sample files located in the test folder
    -----------------------------------------------------------------------------------------------------------------------------------
    --test                        Runs the pipeline as a test run on sample files located on winter server 
    -----------------------------------------------------------------------------------------------------------------------------------
    """.stripIndent()
}
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = "0.7.1"

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}
// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define Configurable Variables
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

params.experiment = false
params.prefix = false
params.sample = false
params.strand = false
params.merge = false
params.input = false
params.unalign = false
params.reference = false
params.annotation = false
params.indexing = false
params.genotype = false
params.output = false
params.email = false
params.name = false
//##TODO(iaguilar): what is params raw? not described in documentation (Doc ######)
//##TODO(iaguilar): raw seems to mean that data requires trimming (Doc ######)
//##TODO(iaguilar): since the pipeline can detect if trimming is necessary, there should not be an option to direct trimming  (Doc ######)
params.raw = false
params.ercc = false
params.fullCov = false
params.test = false
params.small_test = false
params.k_lm = false
params.k_sd = false
//##TODO(iaguilar): Check that no parameter definition is missing (Dev ######)
//  ##TODO(iaguilar): seems that params.scripts is missing (Dev ######)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Validate Inputs
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//##// MANDATORY PARAMS BLOCK
// Sample Selection Validation
if (!params.sample || (params.sample != "single" && params.sample != "paired")) {
    exit 1, "Sample Type Not Provided or Invalid Choice. Please Provide a Valid Sample Type"
}

// Strand Selection Validation
if (!params.strand || (params.strand != "forward" && params.strand != "reverse" && params.strand != "unstranded")) {
    exit 1, "Strand Type Not Provided or Invalid Choice. Please Provide a Valid Strand Type"
}

// Reference Selection Validation
if (!params.reference) {
    exit 1, "Error: enter hg19 or hg38, mm10 for mouse, or rn6 for rat as the reference."
}
if (params.reference == "hg19" || params.reference == "hg38" ) {
    params.reference_type = "human"
}
if (params.reference == "mm10") {
    params.reference_type = "mouse"
}
if (params.reference == "rn6") {
    params.reference_type = "rat"
}

//##// OPTIONAL PARAMS BLOCK
// Annotation Path Validation
if (!params.annotation) {
    params.annotations = "./Annotation"
}

// Indexing Path Validation
if (!params.indexing && !params.test) {
    params.indexing = "${params.annotations}"
}

// Genotype Path Validation
if (!params.genotype) {
    params.genotypes = "./Genotyping"
}

// Experiment/Workflow Name Validation
if (!params.name) {
    if (!params.experiment) {
        workflow.runName = "RNAsp_run"
        params.experiments = "Jlab_experiment"
    }
    if (params.experiment) {
        workflow.runName = params.experiment
        params.experiments = params.experiment
    }
}

// Prefix
if (!params.prefix) {
    params.experiment_prefix = "pref"
}
if (params.prefix) {
    params.experiment_prefix = params.prefix
//  params.experiment_prefix = "Yas"
}

// External Script Path Validation
//##TODO(iaguilar): This param was not defined neither in help nor in the variable definition block (Dev ######)
//##TODO(iaguilar): Comment in original devp was: "It's the directory where the shell files are located at. You only need to specify it if you cloned this repository somewhere else"; since this scripts folder will be versioned with the NF pipeline, there is no need to allow it to be on another directory (Dev ######)
//##TODO(iaguilar): remove this param? only if it is not used by any block, may be replaced by a good old variable to construct paths for script execution (Dev ######)
params.scripts = "./scripts"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Core Options
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//##TODO(iaguilar): Add brief descriptions for what are the cores used (parallelization level, by sample, by chunck, etc) (Doc ######)
params.cores = '1'
params.ercc_cores = '8'
params.trimming_cores = '8'
params.hisat_cores = '8'
params.samtobam_cores = '8'
params.featurecounts_cores = '8'
params.alignments_cores = '8'
params.salmon_cores = '8'
params.counts_cores = '5'
params.coverage_cores = '5'
params.expressedregion_cores = '8'


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Input Path Options
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Base Input Path Options
//##TODO(iaguilar): Explain how the input files must be ordered or located inside the input data directory (Doc ######)
//##TODO(iaguilar): INCLUDE IN THE README: a tree with file strcture description (Doc ######)
if (!params.input && !params.test && !params.small_test) {
    if (params.input) {
        params.inputs = "${params.input}"
    }
    else {
        params.inputs = "./input"
// originally was        params.inputs = "./inputs"
// but it does not seem to be a "inputs" dir in the repo, only "input", without the s
    }
}

// Testing Paths for Rat
//##TODO(iaguilar): Complete with data for rat, by vlarios (Dev ######)
if (params.reference_type == "rat") {
    exit 1, "There are no sample files for rat. Please add sample files and then run again"
}

// Real File Test Paths (Human & Mouse)
//##TODO(iaguilar): This could change to a general test param that runs on every species and genome version (Dev ######)
if (params.test && !params.small_test) {
    if (params.reference_type == "mouse") {
        params.inputs = "/media/genomics/disco3/dataLieber/raw/${params.reference_type}/${params.sample}"
    }
//##TODO(iaguilar): Check what does raw mean (Dev ######)
//##TODO(iaguilar): It means data in need of trimming (Dev ######)
//##TODO(iaguilar): These paths should point to the test data created by Lieber or vlarios (Dev ######)
//##TODO(iaguilar): A minimal test-data directory should be versionated with the pipeline (Dev ######)
    if (params.reference_type == "human") {
        if (params.raw) {
            params.inputs = "/media/genomics/disco3/dataLieber/raw/human/paired"
        }
        if (!params.raw) {
            params.inputs = "/media/genomics/disco3/dataLieber/human"
        }
    }
}

// Dummy File Test Paths
//##TODO(iaguilar): this seems to complex to have to specifiy the test and the type of stranded data (Dev ######)
//##TODO(iaguilar): One or two test runs with multiple species, genome versions, and minimal data should be defined in a single flag (like: "quick_test" and "full_test"  (Dev ######)
if (params.small_test && !params.test) {
    if (params.strand != "unstranded") {
        if (params.merge) {
            params.inputs = "./test/merge/${params.sample}/stranded"
        }
        if (!params.merge) {
            params.inputs = "./test/${params.sample}/stranded"
        }
    }
    if (params.strand == "unstranded") {
        if (params.merge) {
            params.inputs = "./test/merge/${params.sample}/unstranded"
        }
        if (!params.merge) {
            params.inputs = "./test/${params.sample}/unstranded"
        }
    }
}

// Conflicting Test Options
//##TODO(iaguilar): Modify for quick_test vs full_test (Dev ######)
if (params.small_test && params.test) {
    exit 1, "You've selected 'small_test' and 'test' ... Please choose one and run again"
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Strand Option Parameters
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//##TODO(iaguilar): Check if kallisto parameters belong on the stranded option params block, or deserve their own block (Doc ######)
if (!params.k_lm) {
    params.k_length_mean = 200
}
if (!params.k_sd) {
    params.k_standard_deviation = 30
}
if (params.k_lm){
    params.k_length_mean = params.k_lm
}
if (params.k_sd) {
    params.k_standard_deviation = params.k_sd
}
//##TODO(iaguilar): Expand what is every param used for and why is this value assigned (Doc ######)
if (params.strand == "unstranded") {
    if (params.sample == "single") {
        params.kallisto_strand = "--single -l ${params.k_length_mean} -s ${params.k_standard_deviation}"
        params.hisat_strand = ""
        params.feature_strand = "0"
        params.trim_sample = "SE"
    }
    if (params.sample == "paired") {
        params.kallisto_strand = ""
        params.hisat_strand = ""
        params.feature_strand = "0"
        params.trim_sample = "PE"
    }
}
if (params.strand == "forward") {
    if (params.sample == "single") {
        params.kallisto_strand = "--single --fr-stranded -l ${params.k_length_mean} -s ${params.k_standard_deviation}"
        params.hisat_strand = "--rna-strandness F"
        params.feature_strand = "1"
        params.trim_sample = "SE"
    }
    if (params.sample == "paired") {
        params.kallisto_strand = "--fr-stranded"
        params.hisat_strand = "--rna-strandness FR"
        params.feature_strand = "1"
        params.trim_sample = "PE"
    }
}
if (params.strand == "reverse") {
    if (params.sample == "single") {
        params.kallisto_strand = "--single --rf-stranded -l ${params.k_length_mean} -s ${params.k_standard_deviation}"
        params.hisat_strand = "--rna-strandness R"
        params.feature_strand = "2"
        params.trim_sample = "SE"
    }
    if (params.sample == "paired") {
        params.kallisto_strand = "--rf-stranded"
        params.hisat_strand = "--rna-strandness RF"
        params.feature_strand = "2"
        params.trim_sample = "PE"
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Output File Paths (Merging, Paired/Single, Stranded  Combinations)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if (params.output) {
    params.basedir = "${params.output}"
//##TODO(iaguilar): Check if this is necessary or it is best to manage the ./Annotations dir by default (Dev ######)
    params.index_out = "${params.annotations}"
}
if (!params.output) {
//##TODO(iaguilar): What is params.production_baseout used for? (Doc ######)
    params.production_baseout = "."
//##TODO(iaguilar): What is params.test_baseout used for? (Doc ######)
    params.test_baseout = "."
    if (params.test) {
        params.index_out = "${params.test_baseout}/Annotation"
        if (params.merge) {
            params.basedir = "${params.test_baseout}/results/${params.reference_type}/${params.reference}/${params.sample}/merge"
        }
        if (!params.merge) {
            params.basedir = "${params.test_baseout}/results/${params.reference_type}/${params.reference}/${params.sample}"
        }
    }
    if (!params.test) {
        params.index_out = "${params.production_baseout}/Annotation"
        if (params.merge) {
            params.basedir = "${params.production_baseout}/results/${params.reference_type}/${params.reference}/${params.sample}/merge"
        }
        if (!params.merge) {
		params.basedir = "${params.production_baseout}/results/${params.reference_type}/${params.reference}/${params.sample}"
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define External Scripts
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//##TODO(iaguilar): This should live in an external config file... (Dev ######)
infer_strandness = file("${params.scripts}/step3b_infer_strandness.R")
prep_bed = file("${params.scripts}/prep_bed.R")
bed_to_juncs = file("${params.scripts}/bed_to_juncs.py")
fullCov_file = file("${params.scripts}/create_fullCov_object.R")
//TODO(iaguilar) change _file for _script in the variable expressedRegions_file (###Dev)
expressedRegions_file = file("${params.scripts}/step9-find_expressed_regions.R")
check_R_packages_script = file("${params.scripts}/check_R_packages.R")

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define Reference Paths/Scripts + Reference Dependent Parameters
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ERCC
if (params.ercc) {
    erccidx = file("${params.annotations}/ERCC/ERCC92.idx")
}

// ERCC Concentrations
ercc_actual_conc = file("${params.annotations}/ercc_actual_conc.txt")

if (params.reference == "hg38") {
    
    // Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz"
    params.fa_gz = "GRCh38.primary_assembly.genome.fa.gz"
    params.fa = "GRCh38.primary_assembly.genome.fa"
    params.hisat_prefix = "hisat2_GRCh38primary"
    params.hisat_assembly = "GENCODE/GRCh38_hg38/assembly"

    // Step 4: gencode gtf
//##TODO(iaguilar): Check if gtf_link and gtf_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.gencode_gtf_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"
    params.gencode_gtf_gz = "gencode.v25.annotation.gtf.gz"
    params.gencode_gtf = "gencode.v25.annotation.gtf"
    params.feature_output_prefix = "Gencode.v25.hg38"

    // Step 5: python coverage
//##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
    chr_sizes = file("${params.annotations}/chrom_sizes/hg38.chrom.sizes.gencode")

    // Step 6: salmon
//##TODO(iaguilar): Explain why step 6 is enabled if reference is hg38...  (Doc ######)
    params.step6 = true
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.tx_fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.transcripts.fa.gz"
    params.tx_fa_gz = "gencode.v25.transcripts.fa.gz"
    params.tx_fa = "gencode.v25.transcripts.fa"
    params.salmon_prefix = "salmon_0.8.2_index_gencode.v25.transcripts"
    params.salmon_assembly = "GENCODE/GRCh38_hg38/transcripts"

    // Step 7: Make R objects    
    junction_annotation_gencode = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_hg38_gencode_v25.rda")
    junction_annotation_ensembl = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_hg38_ensembl_v85.rda")
    junction_annotation_genes = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_hg38_refseq_grch38.rda")
    feature_to_tx_gencode = Channel.fromPath("${params.annotations}/junction_txdb/feature_to_Tx_hg38_gencode_v25.rda")
    feature_to_tx_ensembl = Channel.fromPath("${params.annotations}/junction_txdb/feature_to_Tx_ensembl_v85.rda")
    create_counts = file("${params.scripts}/create_count_objects-human.R")

    // Step 8: call variants
//##TODO(iaguilar): Explain why step 8 is enabled if reference is hg38...  (Doc ######)
    params.step8 = true
//##TODO(iaguilar): Explain the need to define the channel from this block  (Doc ######)
    Channel
      .fromPath("${params.genotypes}/common_missense_SNVs_hg38.bed")
      .set{ snvbed }

}
if (params.reference == "hg19") {
    
    // Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
    params.fa_gz = "GRCh37.primary_assembly.genome.fa.gz"
    params.fa = "GRCh37.primary_assembly.genome.fa"
    params.hisat_prefix = "hisat2_GRCh37primary"
    params.hisat_assembly = "GENCODE/GRCh37_hg19/assembly"
    
    // Step 4: gencode gtf
//##TODO(iaguilar): Check if gtf_link and gtf_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.gencode_gtf_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz"
    params.gencode_gtf_gz = "gencode.v25lift37.annotation.gtf.gz"
    params.gencode_gtf = "gencode.v25lift37.annotation.gtf"
    params.feature_output_prefix = "Gencode.v25lift37.hg19"

    // Step 5: python coverage
//##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
    chr_sizes = file("${params.annotations}/chrom_sizes/hg19.chrom.sizes.gencode")

    // Step 6: salmon
//##TODO(iaguilar): Explain why step 6 is enabled if reference is hg19...  (Doc ######)
    params.step6 = true
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.tx_fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.transcripts.fa.gz"
    params.tx_fa_gz = "gencode.v25lift37.transcripts.fa.gz"
    params.tx_fa = "gencode.v25lift37.transcripts.fa"
    params.salmon_prefix = "salmon_0.8.2_index_gencode.v25lift37.transcripts"
    params.salmon_assembly = "GENCODE/GRCh37_hg19/transcripts"

    // Step 7: Make R objects
    junction_annotation_gencode = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_hg19_gencode_v25lift37.rda")
    junction_annotation_ensembl = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_hg19_ensembl_v75.rda")
    junction_annotation_genes = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_hg19_refseq_grch37.rda")
    feature_to_tx_gencode = Channel.fromPath("${params.annotations}/junction_txdb/feature_to_Tx_hg19_gencode_v25lift37.rda")
    feature_to_tx_ensembl = Channel.fromPath("${params.annotations}/junction_txdb/feature_to_Tx_ensembl_v75.rda")
    create_counts = file("${params.scripts}/create_count_objects-human.R")

    // Step 8: call variants
//##TODO(iaguilar): Explain why step 8 is enabled if reference is hg19...  (Doc ######)
    params.step8 = true
//##TODO(iaguilar): Explain the need to define the channel from this block  (Doc ######)
    Channel
      .fromPath("${params.genotypes}/common_missense_SNVs_hg19.bed")
      .set{ snvbed }

}
if (params.reference == "mm10") {

    // Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/GRCm38.primary_assembly.genome.fa.gz"
    params.fa_gz = "GRCm38.primary_assembly.genome.fa.gz"
    params.fa = "GRCm38.primary_assembly.genome.fa"
    params.hisat_prefix = "GRCm38_mmhisat2_GRCm38primary"
    params.hisat_assembly = "GENCODE/GRCm38_mm10/assembly"

    // Step 4: gencode gtf
//##TODO(iaguilar): Check if gtf_link and gtf_gz are not redundant since link includes fa_gaz value (Dev ######)
    params.gencode_gtf_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.annotation.gtf.gz"
    params.gencode_gtf_gz = "gencode.vM11.annotation.gtf.gz"
    params.gencode_gtf = "gencode.vM11.annotation.gtf"
    params.feature_output_prefix = "Gencode.M11.mm10"

    // Step 5: python coverage
//##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
    chr_sizes = file("${params.annotations}/chrom_sizes/mm10.chrom.sizes.gencode")
    
    // Step 6: salmon
//##TODO(iaguilar): Explain why step 6 is enabled if reference is mm10...  (Doc ######)
    params.step6 = true
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gz value (Dev ######)
    params.tx_fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.transcripts.fa.gz"
    params.tx_fa_gz = "gencode.vM11.transcripts.fa.gz"
    params.tx_fa = "gencode.vM11.transcripts.fa"
    params.salmon_prefix = "salmon_0.8.2_index_gencode.vM11.transcripts"
    params.salmon_assembly = "GENCODE/GRCm38_mm10/transcripts"

    // Step 7: Make R objects
    junction_annotation_gencode = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_mm10_gencode_vM11.rda")
    junction_annotation_ensembl = Channel.fromPath("${params.annotations}/junction_txdb/junction_annotation_mm10_ensembl_v86.rda")
    create_counts = file("${params.scripts}/create_count_objects-mouse.R")

    // Step 8: call variants
//##TODO(iaguilar): Explain why step 8 is disabled if reference is mm10...  (Doc ######)
    params.step8 = false
}
if (params.reference == "rn6") {

    // Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gz value (Dev ######)
    params.fa_link = "ftp://ftp.ensembl.org/pub/release-86/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
    params.fa_gz = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
    params.fa = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
    params.hisat_prefix = "hisat2_Rnor6.0toplevel"
    params.hisat_assembly = "ensembl/Rnor_6.0"
    
    
    // Step 4: gencode gtf (ensembl for rn6)
//##TODO(iaguilar): Check if gtf_link and gtf_gz are not redundant since link includes gtf_gz value (Dev ######)
    params.gencode_gtf_link = "ftp://ftp.ensembl.org/pub/release-86/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.86.gtf.gz"
    params.gencode_gtf_gz = "Rattus_norvegicus.Rnor_6.0.86.gtf.gz"
    params.gencode_gtf = "Rattus_norvegicus.Rnor_6.0.86.gtf"
    params.feature_output_prefix = "Rnor_6.0.86"

    // Step 5: python coverage
//##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
    chr_sizes = file("${params.annotations}/chrom_sizes/rn6.chrom.sizes.ensembl")
    
    // Step 6: Salmon
//##TODO(iaguilar): Explain why step 6 is enabled if reference is rn6...  (Doc ######)
    params.step6 = false


    //Step 8: call variants
//##TODO(iaguilar): Explain why step 8 is enabled if reference is rn6...  (Doc ######)
    params.step8 = false

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the Prefix Functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_merging_prefix = { file -> file.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_read1_A1s2o2s1A" - "_read2_A1s2o2s1A" }

def get_prefix = { file -> file.name.toString().tokenize('.')[0] }

def get_paired_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_1_A1s2o2s1A" - "_2_A1s2o2s1A" }

def get_summary_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_summary_A1s2o2s1A" }

def get_summary_paired_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_summary_A1s2o2s1A" + "_A1s2o2s1A" - "_1_A1s2o2s1A" - "_2_A1s2o2s1A" }

def get_TR_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_TR_A1s2o2s1A" }

def get_TR_paired_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_TR_A1s2o2s1A" + "_A1s2o2s1A" - "_1_A1s2o2s1A" - "_2_A1s2o2s1A" }

def get_TNR_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_TNR_A1s2o2s1A" } 

def get_TNR_paired_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_TNR_A1s2o2s1A" + "_A1s2o2s1A" - "_1_A1s2o2s1A" - "_2_A1s2o2s1A" }

def get_single_trimmed_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_trimmed_A1s2o2s1A" }

def get_paired_trimmed_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_reverse_paired_A1s2o2s1A" - "_reverse_unpaired_A1s2o2s1A" - "_forward_paired_A1s2o2s1A" - "_forward_unpaired_A1s2o2s1A" + "_A1s2o2s1A" - "_trimmed_A1s2o2s1A" }

def get_hisat_prefix = { file -> file.name.toString().tokenize('.')[0] + "_A1s2o2s1A" - "_hisat_out_A1s2o2s1A" }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Summary of Defined Variables
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Header log info
log.info "============================================================="
log.info " LIBD-RNAseq : Multi-Input RNA-Seq Best Practice v${version}"
log.info "============================================================="
def summary = [:]
summary['Run Name']            = workflow.runName
summary['Reference Type']      = params.reference_type
summary['Sample']              = params.sample
summary['Reference']           = params.reference
summary['Strand']              = params.strand
summary['Annotations']         = params.annotations
summary['Genotypes']           = params.genotypes
summary['Input']               = params.inputs
if(params.ercc) summary['ERCC Index'] = erccidx
if(params.experiments) summary['Experiment'] = params.experiments
if(params.email) summary['E-mail Address'] = params.email
if(params.merge) summary['Merge'] = "True"
if(params.unalign) summary['Align'] = "True"
if(params.fullCov) summary['Full Coverage'] = "True"
if(params.test) summary['Test'] = "True"
summary['Output dir']          = params.basedir
summary['Working dir']         = workflow.workDir
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================="

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN PIPELINE
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/*
 * Step Ia: GENCODE Assembly FA
 */

params.hisat_idx_output = "${params.index_out}/${params.hisat_assembly}"

// get the number of *fa files in the reference directory
Channel
  .fromPath("${params.hisat_idx_output}/fa/*.fa")
  .set{ assembly_fa_trigger }
  
// Use the number of fa files in the reference directory to set the aseembly channel, efectively skipping this step when running the pipeline
// If the count of files in the assembly_fa_trigger is more than zero, the assembly_fa_download channel wont be set
// else assembly_fa_download channel is set, and pipeline starts by downloading the ref genome
assembly_fa_trigger.count().filter{ it == 0 }.set{ assembly_fa_download }

process pullGENCODEassemblyfa {

    echo true
    tag "Downloading Assembly FA File: ${params.fa}"
    publishDir "${params.hisat_idx_output}/fa",'mode':'copy'

	input:
	val(assembly_fa_val) from assembly_fa_download

    output:
    file("${params.fa}") into reference_fa

    script:
    fa_gz_file_link = "${params.fa_link}"
    fa_gz = "${params.fa_gz}"
	fa_file = "${params.fa}"
	fa_file_build = "${fa_file}_build"
	"""
		## Download and unzip the ref genome
		## *_build lock added as intermediate file to be able to quickly find if a process was interrupted in the workdir
		wget "$fa_gz_file_link" \
		&& gunzip -c "$fa_gz" > "$fa_file_build" \
		&& mv "$fa_file_build" "$fa_file"
	"""
}

Channel
  .fromPath("${params.hisat_idx_output}/fa/${params.fa}")
  .mix(reference_fa)
  .toSortedList()
  .flatten()
  .distinct()
  .into{ reference_assembly; variant_assembly }

/*
 * Step Ib: Build HISAT Index
 */

 // Another trigger to check if the next step is necessary
Channel
  .fromPath("${params.hisat_idx_output}/index/${params.hisat_prefix}.*.ht2")
  .set{ hisat_builds_trigger }

  // Check if less than 8 files were found in the previous trigger definition
 // If the count of files in the trigger is less than eigth, the channel wont be set
 // else the channel is set, and pipeline performs this step
hisat_builds_trigger.count().filter{ it < 8 }.set{ hisat_trigger_build }

process buildHISATindex {

    echo true
    tag "Building HISAT2 Index: ${params.hisat_prefix}"
    publishDir "${params.hisat_idx_output}/index",'mode':'copy'

    input:
    val(hisat_trigger_val) from hisat_trigger_build
    file reference_fasta from reference_assembly

    output:
    file("${params.hisat_prefix}.*") into hisat_index_built

    script:
    prefix = "${params.hisat_prefix}"
	// Here, number of threads in the -p option does not seem to be assigned by variable (### Dev)
	// $prefix is "the hisat2_index_base  write ht2 data to files with this dir/basename"
    """
		hisat2-build -p 3 $reference_fasta $prefix
    """
}

// Read the hisat index file from path, and/or mix with the channel output from the buildHISATindex process
// This effectively means that whether the index was created or it was already present, the flow can continue
Channel
  .fromPath("${params.hisat_idx_output}/index/${params.hisat_prefix}.*.ht2")
  .mix(hisat_index_built)
  .toSortedList()
  .flatten()
  .distinct()
  .toSortedList()
  .set{ hisat_index }

/*
 * Step IIa: GENCODE GTF Download
 */

// Define the putput dir for download of the GTF reference
params.gencode_gtf_out = "${params.index_out}/RSeQC/${params.reference}"

 // Another trigger to check if the next step is necessary
 // Will look for a .gtf file in the gtf annotation directory.
 // If it finds it, gencode_gtf_download_trigger will be set; If it does not find it, it will be empty
Channel
  .fromPath("${params.gencode_gtf_out}/gtf/*.gtf")
  .set{ gencode_gtf_download_trigger }

 // If the count of files in the trigger is different from zero, the channel wont be set
 // else the channel is set, and pipeline performs this step
gencode_gtf_download_trigger.count().filter{ it == 0 }.set{ gencode_gtf_trigger_download }

process pullGENCODEgtf {

    echo true
    tag "Downloading GTF File: ${params.gencode_gtf}"
    publishDir "${params.gencode_gtf_out}/gtf",'mode':'copy'

    input:
    val(gencode_gtf_bed_val) from gencode_gtf_trigger_download

    output:
    file("${params.gencode_gtf}") into gencode_gtf_downloaded

    script:
    gencode_gtf_link = "${params.gencode_gtf_link}"
    gencode_gtf_gz = "${params.gencode_gtf_gz}"
	gencode_gtf_file = "${params.gencode_gtf}"
	gencode_gtf_file_build = "${gencode_gtf_file}_build"
    """
    wget "$gencode_gtf_link" \
	&& gunzip -c "$gencode_gtf_gz" > "$gencode_gtf_file_build" \
	&& mv "$gencode_gtf_file_build" "$gencode_gtf_file"
    """
	// May add afterscript to clean the .gz? or the .gtf also if it copies the results BEFORE performing afterscript routines
}

// Read the gencode gtf file from path, and/or mix with the channel output from the pullGENCODEgtf process
// This effectively means that whether the gtf was created or it was already present, the flow can continue
Channel
  .fromPath("${params.gencode_gtf_out}/gtf/${params.gencode_gtf}")
  .mix(gencode_gtf_downloaded)
  .toSortedList()
  .flatten()
  .distinct()
  .into{ gencode_gtf; create_counts_gtf; gencode_feature_gtf }

/*
 * Step IIb: Build Bed File
 */

// Define the putput dir for download of the bed file used by http://rseqc.sourceforge.net/#infer-experiment-py
Channel
  .fromPath("${params.gencode_gtf_out}/bed/*")
  .set{ gencode_gtf_bed_trigger }

// Create a trigger to check if the file exist, similar to previous triggers
gencode_gtf_bed_trigger.count().filter{ it == 0 }.set{ gencode_gtf_trigger_bed }

process buildPrepBED {

    echo true
    tag "Building Bed File: ${params.reference}"
    publishDir "${params.gencode_gtf_out}/bed",'mode':'copy'

    input:
    val(gencode_gtf_bed_val) from gencode_gtf_trigger_bed
    file gencode_gtf from gencode_gtf
    file prep_bed from prep_bed
	file check_R_packages_script from check_R_packages_script

    output:
    file("${name}.bed") into bedfile_built

    shell:
    name = "${params.reference}"
    '''
	## Run the script to check for missing rpackages
	Rscript !{check_R_packages_script} \
	&& Rscript !{prep_bed} -f !{gencode_gtf} -n !{name} \
	// does the resulting bed it need to have permissions 775?
	// currently it is created with 664 (-rw-rw-r--)
    '''
}

// Read the bed file from path, and/or mix with the channel output from the bedfile_built process
// This effectively means that whether the bed was created in this run or it was already present, the flow can continue
Channel
  .fromPath("${params.gencode_gtf_out}/bed/*")
  .mix(bedfile_built)
  .toSortedList()
  .flatten()
  .distinct()
  .set{ bedfile }

// for hg38, hg19, and mm10, step 6 is enabled by params.step6 = true 
// during Define Reference Paths/Scripts + Reference Dependent Parameters
if (params.step6) {

    params.salmon_idx_output = "${params.index_out}/${params.salmon_assembly}"

    /*
     * Step IIIa: GENCODE TX FA Download
     */

	// get the number of *fa files in the reference directory
     Channel
      .fromPath("${params.salmon_idx_output}/fa/*.fa")
      .set{ transcript_fa_download }

	// Use the number of fa files in the reference directory to set the aseembly channel, efectively skipping this step when running the pipeline
	// If the count of files in the transcript_download_trigger is more than zero, the transcript_download channel wont be set
	// else transcript_download channel is set, and pipeline starts by downloading the ref transcriptome
    transcript_fa_download.count().filter{ it == 0 }.set{ transcript_download_trigger }

    process pullGENCODEtranscripts {

        echo true
        tag "Downloading TX FA File: ${params.tx_fa}"
        publishDir "${params.salmon_idx_output}/fa",mode:'copy'

        input:
        val(transcript_trigger_val) from transcript_download_trigger

        output:
        file("${params.tx_fa}") into tx_fa_downloaded
		//TODO (iaguilar) put shell code in same style as for fasta reference download (Dev ###)
        script:
        tx_fa_link = "${params.tx_fa_link}"
        tx_fa_gz = "${params.tx_fa_gz}"
        tx_fa = "${params.tx_fa}"
        """
        wget $tx_fa_link
        gunzip $tx_fa_gz
        """
    }

	// Read the salmon reference from path, and/or mix with the channel output from the pullGENCODEtranscripts process
	// This effectively means that whether the reference was created or it was already present, the flow can continue
     Channel
      .fromPath("${params.salmon_idx_output}/fa/${params.tx_fa}")
      .mix(tx_fa_downloaded)
      .toSortedList()
      .flatten()
      .distinct()
      .set{ transcript_fa }

    /*
     * Step IIIb: Salmon Transcript Build
     */

	// get the number of * files in the reference directory
     Channel
      .fromPath("${params.salmon_idx_output}/salmon/${params.salmon_prefix}/*")
      .set{ salmon_build_trigger }

	// Use the number of files in the reference directory to set the salmon trigger channel
	// If the count of files in the salmon_build_trigger is more than zero, the transcript_download channel wont be set
	// else transcript_download channel is set, and pipeline starts by downloading the ref transcriptome
    salmon_build_trigger.count().filter{ it == 0 }.set{ salmon_trigger_build}

    process buildSALMONindex {

        echo true
        tag "Building Salmon Index: ${params.salmon_prefix}"
        publishDir "${params.salmon_idx_output}/salmon",mode:'copy'

        input:
        file tx_file from transcript_fa
        val(salmon_trigger_val) from salmon_trigger_build

        output:
        file("${params.salmon_prefix}") into salmon_index_built

        script:
        salmon_idx = "${params.salmon_prefix}"
		// in the code execution, -p option is hardcoded to 1 thread
	// TODO (iaguilar): make thread assignation dynamic by using a configurable variable
        """
        salmon index -t $tx_file -i $salmon_idx -p 1 --type quasi -k 31
        """
    }

	// Read the salmon index from path, and/or mix with the channel output from the buildSALMONindex process
	// This effectively means that whether the index was created or it was already present, the flow can continue
     Channel
      .fromPath("${params.salmon_idx_output}/salmon/${params.salmon_prefix}/*")
      .mix(salmon_index_built)
      .toSortedList()
      .flatten()
      .distinct()
      .toSortedList()
      .set{ salmon_index }
}

/*
 * Step A: Run Sample Merging if --merge is specified
 */

// If the --merged flag was used, this block performs fastq.gz file merging
// currently not working, requires testing and debugging
if (params.merge) {

    Channel
      .fromPath("${params.inputs}/*.fastq.gz")
      .toSortedList()
      .flatten()
      .map{file -> tuple(get_merging_prefix(file), file) }
      .groupTuple()
      .ifEmpty{ error "Could not find file pairs for merging"}
      .set{ unmerged_pairs }

      process sampleMerging {

        echo true
        tag "prefix: $merging_prefix | Sample Pair: [ $unmerged_pair ]"
        publishDir "${params.basedir}/merged_fastq",'mode':'copy'

        input:
        set val(merging_prefix), file(unmerged_pair) from unmerged_pairs

        output:
        file "*.fastq.gz" into ercc_merged_inputs, fastqc_merged_inputs
		
		shell:
        '''
		## Use of a local prefiz variable helps to avoid generating dump data in input dir due to fullpaths being captured by merging_prefix NF variable
		local_prefix=`echo !{merging_prefix} | rev | cut -d "/" -f1 | rev`
		## To follow illumina naming conventions,_read1 and _read2 will be changed to R1 and R2
		read1="${local_prefix}_read1.fastq.gz"
		read2="${local_prefix}_read2.fastq.gz"
        zcat "${read1}" "${read2}" \
		| gzip -c > "${local_prefix}.fastq.gz"
        '''
    }
}

/*
 * Step B: Run the ERCC process if the --ercc flag is specified
 */

if (params.ercc) {

    if (params.merge) {

        if (params.sample == "single") {

            ercc_merged_inputs
              .flatten()
              .toSortedList()
              .flatten()
              .map{ file -> tuple(get_prefix(file), file) }
              .ifEmpty{ error "Could not find Channel for Merged Single Sample Files for ERCC" }
              .set{ ercc_inputs }
        }
        if (params.sample == "paired") {

            ercc_merged_inputs
              .flatten()
              .toSortedList()
              .flatten()
              .map{file -> tuple(get_paired_prefix(file), file) }
              .groupTuple()
              .ifEmpty{ error "Could not find Channel for Merged Paired Sample Files for ERCC" }
              .set{ ercc_inputs }
        }
    }
    if (!params.merge) {

        if (params.sample == "single") {

            Channel
              .fromPath("${params.inputs}/*.fastq.gz")
              .flatten()
              .toSortedList()
              .flatten()
              .map{ file -> tuple(get_prefix(file), file) }
              .ifEmpty{ error "Could not Find Unmerged Sample Files for ERCC"}
              .set{ ercc_inputs }
        }
        if (params.sample == "paired") {

            Channel
              .fromPath("${params.inputs}/*.fastq.gz")
              .flatten()
              .toSortedList()
              .flatten()
              .map{ file -> tuple(get_paired_prefix(file), file) }
              .groupTuple()
              .ifEmpty{ error "Could not Find Unmerged Sample Files for ERCC"}
              .set{ ercc_inputs }
        }
    }

     process sampleERCC {

        echo true
        tag "Prefix: $ercc_prefix | Sample: [ $ercc_input ]"
        publishDir "${params.basedir}/ercc/${ercc_prefix}",'mode':'copy'

        input:
        file erccidx from erccidx
        set val(ercc_prefix), file(ercc_input) from ercc_inputs

        output:
//        file "*_abundance.tsv" into ercc_abundances
        file("${ercc_prefix}_abundance.tsv") into ercc_abundances

        script:
        ercc_cores = "${params.ercc_cores}"
        strand_option = "${params.kallisto_strand}"
        """
        kallisto quant -i $erccidx -t $ercc_cores -o . $strand_option $ercc_input \
		&& cp abundance.tsv ${ercc_prefix}_abundance.tsv
        """
    }
}

if (params.merge) {

    if (params.sample == "single") {

        fastqc_merged_inputs
          .flatten()
          .toSortedList()
          .flatten()
          .map{file -> tuple(get_prefix(file), file) }
          .ifEmpty{ error "Could not find Channel for Merged Sample Files for FastQC" }
          .into{ fastqc_untrimmed_inputs; adaptive_trimming_fastqs; manifest_creation; salmon_inputs }
    }
    if (params.sample == "paired") {

        fastqc_merged_inputs
          .flatten()
          .toSortedList()
          .flatten()
          .map{file -> tuple(get_paired_prefix(file), file) }
          .groupTuple()
          .ifEmpty{ error "Could not find Channel for Merged Sample Files for FastQC" }
          .into{ fastqc_untrimmed_inputs; adaptive_trimming_fastqs; manifest_creation; salmon_inputs }
    }
}
if (!params.merge) {

    if (params.sample == "single") {

        Channel
          .fromPath("${params.inputs}/*.fastq.gz")
          .flatten()
          .toSortedList()
          .flatten()
          .map{file -> tuple(get_prefix(file), file) }
          .ifEmpty{ error "Could not Find Unmerged Untrimmed Single Sample Files for FastQC"}
          .into{ fastqc_untrimmed_inputs; adaptive_trimming_fastqs; manifest_creation; salmon_inputs }
    }
    if (params.sample == "paired") {

        Channel
          .fromPath("${params.inputs}/*.fastq.gz")
          .flatten()
          .toSortedList()
          .flatten()
          .map{file -> tuple(get_paired_prefix(file), file) }
          .groupTuple()
          .ifEmpty{ error "Could not Find Unmerged Untrimmed Paired Sample Files for FastQC"}
          .into{ fastqc_untrimmed_inputs; adaptive_trimming_fastqs; manifest_creation; salmon_inputs }
    }
}

/*
 * Step C1: Individual Sample Manifest
 */

process sampleIndividualManifest {

    echo true
    tag "Individual Manifest: $manifest_samples $samples_prefix > samples.manifest.${samples_prefix}"
    publishDir "${params.basedir}/manifest",'mode':'copy'

    input:
    set val(samples_prefix), file(manifest_samples) from manifest_creation

    output:
    file "samples.manifest.${samples_prefix}" into individual_manifests

    script:
    """
    printf "${manifest_samples} ${samples_prefix}\n" >> "samples.manifest.${samples_prefix}"
    """
}

individual_manifests
  .flatten()
  .collect()
  .set{ individual_manifest_files }

/*
 * Step C2: Sample Manifest
 */

process sampleManifest {

    echo true
    tag "Aggregate Manifest: $individual_manifests > samples.manifest"
    publishDir "${params.basedir}/manifest",mode:'copy'

    input:
    file individual_manifests from individual_manifest_files

    output:
    file "samples.manifest" into counts_samples_manifest, fullCov_samples_manifest

    script:
    """
    cat ${individual_manifests} > "samples.manifest"
    """
}

/*
 * Step 1: Untrimmed Quality Report
 */

process sampleQualityUntrimmed {

    echo true
    tag "Prefix: $untrimmed_prefix | Sample: [ $fastqc_untrimmed_input ]"
    publishDir "${params.basedir}/FastQC/Untrimmed",mode:'copy'

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
    }
    if (params.sample == "paired") {
        copy_command = "cp ${untrimmed_prefix}_1_fastqc/summary.txt ${untrimmed_prefix}_1_summary.txt && cp ${untrimmed_prefix}_2_fastqc/summary.txt ${untrimmed_prefix}_2_summary.txt"
        data_command = "cp ${untrimmed_prefix}_1_fastqc/fastqc_data.txt ${untrimmed_prefix}_1_fastqc_data.txt && cp ${untrimmed_prefix}_2_fastqc/fastqc_data.txt ${untrimmed_prefix}_2_fastqc_data.txt"
    }
    """
    fastqc $fastqc_untrimmed_input --extract
    $copy_command
    $data_command
    """
}

/*
 * Step 2a: Adaptive Trimming
 * Merge the quality reports channel and trimming input channels
 */

if (params.sample == "single") {

    quality_reports
      .flatten()
      .map{ file -> tuple(get_summary_prefix(file), file) }
      .join(adaptive_trimming_fastqs)
      .ifEmpty{ error "Cannot Find Combined Quality and Trimming Channel for Single Adaptive Trimming" }
      .set{ adaptive_trimming_single_inputs }

    process sampleAdaptiveTrimSingleReads {

      echo true
      tag "Prefix: $single_adaptive_prefix : Sample: [ $single_adaptive_fastq | $single_adaptive_summary ]"
      publishDir "${params.basedir}/Adaptive_Trim",'mode':'copy'

      input:
      set val(single_adaptive_prefix), file(single_adaptive_summary), file(single_adaptive_fastq) from adaptive_trimming_single_inputs

      output:
      file "*" into trimming_fastqs, no_trimming_fastqs

      shell:
      single_quality_report = single_adaptive_prefix.toString() + "_summary.txt"
      single_trimming_input = single_adaptive_prefix.toString() + ".fastq.gz"
      '''
      export result=$(grep "Adapter Content" !{single_quality_report} | cut -f1)
      if [ $result == "FAIL" ] ; then
          mv !{single_trimming_input} "!{single_adaptive_prefix}_TR.fastq.gz"
      else
          mv !{single_trimming_input} "!{single_adaptive_prefix}_TNR.fastq.gz"
      fi
      '''
    }
}
/* if (params.sample == "paired") {

    quality_reports
      .flatten()
      .map{ file -> tuple(get_summary_paired_prefix(file), file) }
      .groupTuple()
      .join(adaptive_trimming_fastqs)
      .ifEmpty{ error "Cannot Find Combined Quality and Trimming Channel for Paired Adaptive Trimming" }
      .set{ adaptive_trimming_paired_inputs }

    process sampleAdaptiveTrimPairedReads {

      echo true
      tag "Prefix: $paired_adaptive_prefix | Sample: [ $paired_adaptive_fastq | $paired_adaptive_summary ]"
      publishDir "${params.basedir}/Adaptive_Trim",mode:'copy'

      input:
      set val(paired_adaptive_prefix), file(paired_adaptive_summary), file(paired_adaptive_fastq) from adaptive_trimming_paired_inputs

      output:
      file "*" into trimming_fastqs, no_trimming_fastqs

      shell:
      quality_report_1 = paired_adaptive_prefix.toString() + "_1_summary.txt"
      quality_report_2 = paired_adaptive_prefix.toString() + "_2_summary.txt"
      trimming_input_1 = paired_adaptive_prefix.toString() + "_1.fastq.gz"
      trimming_input_2 = paired_adaptive_prefix.toString() + "_2.fastq.gz"
      adaptive_out_prefix_1 = paired_adaptive_prefix.toString() + "_1"
      adaptive_out_prefix_2 = paired_adaptive_prefix.toString() + "_2"

      '''
      export result1=$(grep "Adapter Content" !{quality_report_1} | cut -c1-4)
      export result2=$(grep "Adapter Content" !{quality_report_2} | cut -c1-4)
      if [ $result1 == "FAIL" || $result2 == "FAIL"] ; then
          cp !{trimming_input_1} "!{adaptive_out_prefix_1}_TR.fastq.gz"
          cp !{trimming_input_2} "!{adaptive_out_prefix_2}_TR.fastq.gz"
      else
          cp !{trimming_input_1} "!{adaptive_out_prefix_1}_TNR.fastq.gz"
          cp !{trimming_input_2} "!{adaptive_out_prefix_2}_TNR.fastq.gz"
      fi
      '''
    }
} */

/*
 * Modify the Trimming Input Channel 
 */

// if (params.sample == "single") {

    // trimming_fastqs
      // .flatten()
      // .filter{ file -> file.name.toString() =~ /_TR.*/ }
      // .map{ file -> tuple(get_TR_prefix(file), file) }
      // .set{ trimming_inputs }

    // no_trimming_fastqs
      // .flatten()
      // .filter{ file -> file.name.toString() =~ /_TNR.*/ }
      // .map{ file -> tuple(get_TNR_prefix(file), file) }
      // .set{ no_trim_fastqs }
  // }

  // if (params.sample == "paired") {

    // trimming_fastqs
      // .flatten()
      // .filter{ file -> file.name.toString() =~ /_TR.*/ }
      // .toSortedList()
      // .flatten()
      // .map{ file -> tuple(get_TR_paired_prefix(file), file) }
      // .groupTuple()
      // .set{ trimming_inputs }

    // no_trimming_fastqs
      // .flatten()
      // .filter{ file -> file.name.toString() =~ /_TNR.*/ }
      // .toSortedList()
      // .flatten()
      // .map{ file -> tuple(get_TNR_paired_prefix(file), file) }
      // .groupTuple()
      // .set{ no_trim_fastqs }
// }

/*
 * Step 2b: Trimming 
 */

/* process sampleTrimming {

    echo true
    tag "Prefix: $trimming_prefix | Sample: [ $trimming_input ]"
    publishDir "${params.basedir}/trimmed_fq",mode:'copy'

    input:
    set val(trimming_prefix), file(trimming_input) from trimming_inputs

    output:
    file "*.fastq.gz" into trimmed_fastqc_inputs, trimmed_hisat_inputs

    script:
    trimming_cores = "${params.trimming_cores}"
    sample_option = "${params.trim_sample}"
    if (params.sample == "single") {
        output_option = "${trimming_prefix}_trimmed.fastq.gz"
    }
    if (params.sample == "paired") {
        output_option = "${trimming_prefix}_trimmed_forward_paired.fastq.gz ${trimming_prefix}_trimmed_forward_unpaired.fastq.gz ${trimming_prefix}_trimmed_reverse_paired.fastq.gz ${trimming_prefix}_trimmed_reverse_unpaired.fastq.gz"
    }
    """
    java -Xmx512M \
    -jar /usr/local/bin/trimmomatic-0.36.jar \
    $sample_option \
    -threads $trimming_cores \
    -phred33 \
    $trimming_input \
    $output_option \
    ILLUMINACLIP:/usr/local/TruSeq2-PE.fa:2:30:10:1 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:75
    """
} */

/*
 * Step 2c: Run FastQC Quality Check on Trimmed Files
 */

/* process sampleQualityTrimmed {
  
    echo true
    tag "$fastqc_trimmed_input"
    publishDir "${params.basedir}/FastQC/Trimmed",mode:'copy'

    input:
    file fastqc_trimmed_input from trimmed_fastqc_inputs

    output:
    file "*"

    script:
    """
    fastqc $fastqc_trimmed_input --extract
    """
} */

/*
 * Step 3a: Hisat Sam File
 */

/* if (params.sample == "single") {

    trimmed_hisat_inputs
      .flatten()
      .map{ file -> tuple(get_single_trimmed_prefix(file), file) }
      .mix(no_trim_fastqs)
      .ifEmpty{ error "Single End Channel for HISAT is empty" }
      .set{ single_hisat_inputs } */


/*     process sampleSingleEndHISAT {

      echo true
      tag "Prefix: $single_hisat_prefix | Sample: $single_hisat_input"
      publishDir "${params.basedir}/HISAT2_out",mode:'copy'

      input:
      file hisat_index from hisat_index
      set val(single_hisat_prefix), file(single_hisat_input) from single_hisat_inputs

      output:
      file "*_hisat_out.sam" into hisat_single_output
      file "*"
      file "*_align_summary.txt" into alignment_summaries

      shell:
      hisat_prefix = "${params.hisat_prefix}"
      strand = "${params.hisat_strand}"
      hisat_cores = "${params.hisat_cores}"
      '''
      hisat2 \
      -p !{hisat_cores} \
      -x !{hisat_prefix} \
      -U !{single_hisat_input} \
      -S !{single_hisat_prefix}_hisat_out.sam !{strand} --phred33 \
      2> !{single_hisat_prefix}_align_summary.txt
      '''
    } */

/*     hisat_single_output
      .flatten()
      .map{ file -> tuple(get_hisat_prefix(file), file) }
      .set{ sam_to_bam_inputs }
} */


/* if (params.sample == "paired") {

    no_trim_fastqs
      .set{ notrim_paired_hisat_inputs }

    process samplePairedEndNoTrimHISAT {

      echo true
      tag "Prefix: $paired_notrim_hisat_prefix | Sample: [ $paired_no_trim_hisat ]"
      publishDir "${params.basedir}/HISAT2_out",mode:'copy'

      input:
      file hisatidx from hisat_index
      set val(paired_notrim_hisat_prefix), file(paired_no_trim_hisat) from notrim_paired_hisat_inputs

      output:
      file "*_hisat_out.sam" into hisat_paired_notrim_output
      file "*"
      file "*_align_summary.txt" into paired_notrim_alignment_summaries

      shell:
      hisat_prefix = "${params.hisat_prefix}"
      strand = "${params.hisat_strand}"
      hisat_cores = "${params.hisat_cores}"
      sample_1_hisat = paired_notrim_hisat_prefix.toString() + "_1_TNR.fastq.gz"
      sample_2_hisat = paired_notrim_hisat_prefix.toString() + "_2_TNR.fastq.gz"
      if (params.unalign) {
          unaligned_opt = "--un-conc ${paired_notrim_hisat_prefix}.fastq"
      }
      if (!params.unalign) {
          unaligned_opt = ""
      }
      '''
      hisat2 \
      -p !{hisat_cores} \
      -x !{hisat_prefix} \
      -1 !{sample_1_hisat} \
      -2 !{sample_2_hisat} \
      -S !{paired_notrim_hisat_prefix}_hisat_out.sam !{strand} --phred33 \
      !{unaligned_opt} \
      2> !{paired_notrim_hisat_prefix}_align_summary.txt
      '''
    } */

/*     trimmed_hisat_inputs
      .flatten()
      .map{ file -> tuple(get_paired_trimmed_prefix(file), file) }
      .groupTuple()
      .set{ trim_paired_hisat_inputs } */


/*     process samplePairedEndTrimmedHISAT {

      echo true
      tag "Prefix: $paired_trimmed_prefix | Sample: $paired_trimmed_fastqs"
      publishDir "${params.basedir}/HISAT2_out",mode:'copy'

      input:
      file hisatidx from hisat_index
      set val(paired_trimmed_prefix), file(paired_trimmed_fastqs) from trim_paired_hisat_inputs

      output:
      file "*_hisat_out.sam" into hisat_paired_trim_output
      file "*"
      file "*_align_summary.txt" into paired_trim_alignment_summaries

      shell:
      hisat_prefix = "${params.hisat_prefix}"
      strand = "${params.hisat_strand}"
      hisat_cores = "${params.hisat_cores}"
      forward_paired = paired_trimmed_prefix.toString() + "_trimmed_forward_paired.fastq.gz"
      reverse_paired = paired_trimmed_prefix.toString() + "_trimmed_forward_paired.fastq"
      forward_unpaired = paired_trimmed_prefix.toString() + "_trimmed_forward_unpaired.fastq.gz"
      reverse_unpaired = paired_trimmed_prefix.toString() + "_trimmed_forward_unpaired.fastq.gz"
      if (params.unalign) {
          unaligned_opt = "--un-conc ${paired_trimmed_prefix}.fastq"
      }
      if (!params.unalign) {
          unaligned_opt = ""
      }
      '''
      hisat2 \
      -p !{hisat_cores} \
      -x !{hisat_prefix} \
      -1 !{forward_paired} \
      -2 !{reverse_paired} \
      -U !{forward_unpaired} , !{reverse_unpaired} \
      -S !{paired_trimmed_prefix}_hisat_out.sam !{strand} --phred33 \
      !{unaligned_opt} \
      2> !${paired_trimmed_prefix}_align_summary.txt
      '''
    } */

/*     hisat_paired_notrim_output
      .mix(hisat_paired_trim_output)
      .flatten()
      .map{ file -> tuple(get_hisat_prefix(file), file) }
      .set{ sam_to_bam_inputs }

    paired_trim_alignment_summaries
      .mix(paired_notrim_alignment_summaries)
      .flatten()
      .set{ alignment_summaries }
} */

/*
 * Step 3b: Sam to Bam 
 */


/* process sampleSamtoBam {

    echo true
    tag "Prefix: $sam_to_bam_prefix | Sample: $sam_to_bam_input"
    publishDir "${params.basedir}/HISAT2_out/sam_to_bam",mode:'copy'

    input:
    set val(sam_to_bam_prefix), file(sam_to_bam_input) from sam_to_bam_inputs

    output:
    set val("${sam_to_bam_prefix}"), file("${sam_to_bam_prefix}*.sorted.bam"), file("${sam_to_bam_prefix}*.sorted.bam.bai") into infer_experiment_inputs, feature_bam_inputs, alignment_bam_inputs, coverage_bam_inputs, full_coverage_bams, count_objects_bam_files, variant_calls_bam

    script:
    original_bam = "${sam_to_bam_prefix}_accepted_hits.bam"
    sorted_bam = "${sam_to_bam_prefix}_accepted_hits.sorted"
    samtobam_cores = "${params.samtobam_cores}"
    """
    samtools view -bh -F 4 $sam_to_bam_input > $original_bam
    samtools sort -@ $samtobam_cores $original_bam -o ${sorted_bam}.bam
    samtools index ${sorted_bam}.bam
    """
} */

/* infer_experiment_inputs
  .combine(bedfile)
  .set{ infer_experiments } */

/*
 * Step 3c: Infer Experiment
 */

/* process sampleInferExperiment {

    echo true
    tag "Prefix: $infer_prefix | Sample: $bam_file | Index: $bam_index"
    publishDir "${params.basedir}/HISAT2_out/infer_experiment",mode:'copy'

    input:
    set val(infer_prefix), file(bam_file), file(bam_index), file(bed_file) from infer_experiments

    output:
    file "*"
    file "*_experiment.txt" into infer_experiment_outputs

    shell:
    '''
    python /usr/local/bin/infer_experiment.py \
    -i !{bam_file} \
    -r !{bed_file} \
    1> !{infer_prefix}_experiment.txt \
    2> !{infer_prefix}_experiment_summary_out.txt
    '''
} */

/* infer_experiment_outputs
  .collect()
  .flatten()
  .toSortedList()
  .set{ infer_experiment_output } */

/*
 * Step 3d: Infer Strandness
 */

/* process sampleInferStrandness {

    echo true
    tag "Sample: $infer_experiment_files"
    publishDir "${params.basedir}/HISAT2_out/infer_strandness/",mode:'copy'

    input:
    file infer_strandness from infer_strandness
    file infer_experiment_files from infer_experiment_output

    output:
    file "*"
    file "inferred_strandness_pattern.txt" into inferred_strand_coverage, inferred_strand_mean_coverage, inferred_strand_objects

    shell:
    inferred_strandness_pattern = "inferred_strandness_pattern.txt"
    '''
    Rscript !{infer_strandness} -p !{inferred_strandness_pattern}
    '''
} */

/* feature_bam_inputs
  .combine(gencode_feature_gtf)
  .set{ feature_counts_inputs } */

/*
 * Step 4a: Feature Counts
 */

/* process sampleFeatureCounts {

    echo true
    tag "Prefix: $feature_prefix | Sample: $feature_bam | Sample Index: $feature_index"
    publishDir "${params.basedir}/Counts",mode:'copy'

    input:
    set val(feature_prefix), file(feature_bam), file(feature_index), file(gencode_gtf_feature) from feature_counts_inputs

    output:
    file "*"
    file "*.counts" into sample_counts

    script:
    if (params.sample == "single") {
        sample_option = ""
    }
    if (params.sample == "paired") {
        sample_option = "-p"
    }
    feature_out = "${feature_prefix}_${params.feature_output_prefix}"
    featurecounts_cores = "${params.featurecounts_cores}"
    feature_strand = "${params.feature_strand}"
    """ 
    featureCounts \
    -s $feature_strand \
    $sample_option \
    -T $featurecounts_cores \
    -a $gencode_gtf_feature \
    -o ${feature_out}_Genes.counts \
    $feature_bam

    featureCounts \
    -s $feature_strand \
    $sample_option \
    -O \
    -f \
    -T $featurecounts_cores \
    -a $gencode_gtf_feature \
    -o ${feature_out}_Exons.counts \
    $feature_bam
    """
} */

/*
 * Step 4b: Primary Alignments
 */

/* process samplePrimaryAlignments {

    echo true
    tag "Prefix: $alignment_prefix | Sample: [ $alignment_bam ]"
    publishDir "${params.basedir}/Counts/junction/primary_aligments",mode:'copy'

    input:
    set val(alignment_prefix), file(alignment_bam), file(alignment_index) from alignment_bam_inputs

    output:
    set val("${alignment_prefix}"), file("${alignment_prefix}.bam"), file("${alignment_prefix}.bam.bai") into primary_alignments

    script:
    alignments_cores = "${params.alignments_cores}"
    """
    samtools view -@ $alignments_cores -bh -F 0x100 $alignment_bam > ${alignment_prefix}.bam
    samtools index ${alignment_prefix}.bam
    """
} */

/*
 * Step 4c: Junctions
 */

/* process sampleJunctions {

    echo true
    tag "Prefix: $junction_prefix | Sample: [ $alignment_bam ]"
    publishDir "${params.basedir}/Counts/junction",mode:'copy'

    input:
    file bed_to_juncs from bed_to_juncs
    set val(junction_prefix), file(alignment_bam), file(alignment_index) from primary_alignments

    output:
    file "*"
    file("*.count")

    shell:
    outjxn = "${junction_prefix}_junctions_primaryOnly_regtools.bed"
    outcount = "${junction_prefix}_junctions_primaryOnly_regtools.count"
    '''
    regtools junctions extract -i 9 -o !{outjxn} !{alignment_bam}
    python !{bed_to_juncs} < !{outjxn} > !{outcount}
    '''
} */

/*
 * Step 5a: Coverage
 */

/* process sampleCoverage {

    echo true
    tag "Prefix: $coverage_prefix | Infer: $inferred_strand | Sample: $sorted_coverage_bam ]"
    publishDir "${params.basedir}/Coverage/wigs",mode:'copy'

    input:
    file inferred_strand from inferred_strand_coverage
    set val(coverage_prefix), file(sorted_coverage_bam), file(sorted_bam_index) from coverage_bam_inputs
    file chr_sizes from chr_sizes

    output:
    set val("${coverage_prefix}"), file("${coverage_prefix}.wig") into wig_files

    shell:
    '''
    export coverage_strand_rule=$(cat !{inferred_strand})
    if [ $coverage_strand_rule == "none" ] ; then
        python /usr/local/bin/bam2wig.py -s !{chr_sizes} -i !{sorted_coverage_bam} -t 4000000000 -o !{coverage_prefix}
    elif [ $coverage_strand_rule == "1++,1--,2+-,2-+" ] ; then
        python /usr/local/bin/bam2wig.py -s !{chr_sizes} -i !{sorted_coverage_bam} -t 4000000000 -o !{coverage_prefix} -d "1++,1--,2+-,2-+"
    elif [ $coverage_strand_rule == "1+-,1-+,2++,2--" ] ; then
        python /usr/local/bin/bam2wig.py -s !{chr_sizes} -i !{sorted_coverage_bam} -t 4000000000 -o !{coverage_prefix} -d "1+-,1-+,2++,2--"
    elif [ $coverage_strand_rule == "++,--" ] ; then
      python /usr/local/bin/bam2wig.py -s !{chr_sizes} -i !{sorted_coverage_bam} -t 4000000000 -o !{coverage_prefix} -d "++,--"
    elif [ $coverage_strand_rule == "+-,-+" ] ; then
        python /usr/local/bin/bam2wig.py -s !{chr_sizes} -i !{sorted_coverage_bam} -t 4000000000 -o !{coverage_prefix} -d "+-,-+"
    fi
    '''
} */

/*
 * Step 5b: Wig to BigWig
 */

/* process sampleWigToBigWig {

    echo true
    tag "Prefix: $wig_prefix | Sample: [ $wig_file ]"
    publishDir "${params.basedir}/Coverage/BigWigs",mode:'copy'

    input:
    set val(wig_prefix), file(wig_file) from wig_files
    file chr_sizes from chr_sizes

    output:
    file "*.bw" into coverage_bigwigs

    shell:
    '''
    /usr/local/bin/wigToBigWig !{wig_file} !{chr_sizes} !{wig_prefix}.bw
    '''

} */

/* coverage_bigwigs
  .collect()
  .flatten()
  .toSortedList()
  .into{ mean_coverage_bigwigs;full_coverage_bigwigs } */

/*
 * Step 5c: Mean Coverage
 */

/* process sampleMeanCoverage {

    echo true
    tag "Samples: [ $mean_coverage_bigwig ]"
    publishDir "${params.basedir}/Coverage/mean",mode:'copy'

    input:
    file inferred_strand_file from inferred_strand_mean_coverage
    file mean_coverage_bigwig from mean_coverage_bigwigs
    file chr_sizes from chr_sizes

    output:
    file "mean*.bw" into mean_bigwigs, expressed_regions_mean_bigwigs

    shell:
    '''
    export coverage_strand_rule=$(cat !{inferred_strand_file})
    if [ $coverage_strand_rule == "none" ] ; then
        /usr/local/bin/wiggletools write mean.wig mean !{mean_coverage_bigwig}
        /usr/local/bin/wigToBigWig mean.wig !{chr_sizes} mean.bw
    else
        /usr/local/bin/wiggletools write mean.forward.wig mean !{mean_coverage_bigwig}
        /usr/local/bin/wigToBigWig mean.forward.wig !{chr_sizes} mean.forward.bw
        /usr/local/bin/wiggletools write mean.reverse.wig mean !{mean_coverage_bigwig}
        /usr/local/bin/wigToBigWig mean.reverse.wig !{chr_sizes} mean.reverse.bw
    fi
    '''
} */

// for hg38, hg19, and mm10, step 6 is enabled by params.step6 = true during Define Reference Paths/Scripts + Reference Dependent Parameters
/* if (params.step6) { */
    /*
     * Step 6: txQuant
     */

/*     process sampleTXQuant {

        echo true
        tag "Prefix: $salmon_input_prefix | Sample: [ $salmon_inputs ]"
        publishDir "${params.basedir}/Salmon_tx/${salmon_input_prefix}",mode:'copy'

        input:
        file salmon_index from salmon_index
        set val(salmon_input_prefix), file(salmon_inputs) from salmon_inputs

        output:
        file "${salmon_input_prefix}/*"
        file("${salmon_input_prefix}_quant.sf") into salmon_quants

        shell:
        salmon_cores = "${params.salmon_cores}"
        salmon_index_prefix = "${params.salmon_prefix}"
        if (params.sample == "single") {
            sample_command = "-r ${salmon_input_prefix}.fastq.gz"
            if (params.strand == "unstranded" ) {
                salmon_strand = "U"
            }
            if (params.strand == "forward" ) {
                salmon_strand = "SF"
            }
            if (params.strand == "reverse" )
                salmon_strand = "SR"
            }
        if (params.sample == "paired") {
            sample_command = "-1 ${salmon_input_prefix}_1.fastq.gz -2 ${salmon_input_prefix}_2.fastq.gz"
            if (params.strand == "unstranded" ) {
                salmon_strand = "IU"
            }
            if (params.strand == "forward" ) {
                salmon_strand = "ISF"
            }
            if (params.strand == "reverse" ) {
                salmon_strand = "ISR"
            }
        } */
/*         '''
        mkdir !{salmon_index_prefix}
        cp !{salmon_index} !{salmon_index_prefix}/.
        salmon quant \
        -i !{salmon_index_prefix} \
        -p !{salmon_cores} \
        -l !{salmon_strand} \
        !{sample_command} \
        -o !{salmon_input_prefix}
        cp !{salmon_input_prefix}/quant.sf !{salmon_input_prefix}_quant.sf
        '''
    }
} */

/*
 * Construct the Counts Objects Input Channel
 */


/* count_objects_bam_files
  .flatten()
  .buffer(size:2,skip:1)
  .flatten()
  .mix(count_objects_quality_reports)
  .mix(count_objects_quality_metrics)
  .mix(alignment_summaries)
  .mix(create_counts_gtf)
  .mix(sample_counts)
  .collect()
  .flatten()
  .toSortedList()
  .set{ counts_objects_channel } */

/* if (params.reference_type == "human" || params.reference_type == "mouse") {

    counts_objects_channel
      .mix(salmon_quants)
      .collect()
      .flatten()
      .toSortedList()
      .set{counts_objects_channel_1}
}
if (params.reference_type == "rat") {

    counts_objects_channel
      .set{counts_objects_channel_1}
}
if (params.ercc) {
        
    counts_objects_channel_1
      .mix(ercc_abundances)
      .collect()
      .flatten()
      .toSortedList()
      .set{ counts_inputs }
}
if (!params.ercc) {

    counts_objects_channel_1
      .collect()
      .flatten()
      .toSortedList()
      .set{ counts_inputs }
} */


/*
 * Construct the Annotation Input Channel
 */

/* junction_annotation_ensembl
  .collect()
  .flatten()
  .toSortedList()
  .set{rat_annotations} */

/* if (params.reference_type == "rat") {
    rat_annotations
      .collect()
      .flatten()
      .toSortedList()
      .set{counts_annotations}
} */

/* if(params.reference_type != "rat") {

    rat_annotations
      .mix(junction_annotation_gencode)
      .collect()
      .flatten()
      .toSortedList()
      .set{mouse_annotations}
}

if(params.reference_type == "mouse") {
    mouse_annotations
      .collect()
      .flatten()
      .toSortedList()
      .set{counts_annotations}
}

if(params.reference_type == "human") {
    mouse_annotations
      .mix(junction_annotation_genes)
      .mix(feature_to_tx_gencode)
      .mix(feature_to_tx_ensembl)
      .collect()
      .flatten()
      .toSortedList()
      .set{counts_annotations}
} */

/*
 * Step 7a: Create Count Objects
 */

/* process sampleCreateCountObjects {

    echo true
    tag "Creating Counts Objects: [ $counts_input ] | Annotations: [ $counts_annotation ]"
    publishDir "${params.basedir}/Count_Objects",mode:'copy'

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
    }
    if (!params.ercc) {
        ercc_bool = "FALSE"
    }
    if (params.sample == "paired") {
        counts_pe = "TRUE"
    }
    if (params.sample == "single") {
        counts_pe = "FALSE"
    }
    if (params.strand == "unstranded") {
        counts_strand = "-s FALSE"
    }
    if (params.strand == "forward") {
        counts_strand = "-s forward"
    }
    if (params.strand == "reverse") {
        counts_strand = "-s reverse"
    }
    counts_cores = "${params.counts_cores}"
    counts_reference = "${params.reference}"
    counts_experiment = "${params.experiments}"
    counts_prefix = "${params.experiment_prefix}"
    counts_dir = "./"
    '''
    Rscript !{create_counts} -o !{counts_reference} -m !{counts_dir} -e !{counts_experiment} -p !{counts_prefix} -l !{counts_pe} -c !{ercc_bool} -t !{counts_cores} !{counts_strand}
    '''
} */

/* if (params.fullCov) {

    full_coverage_bams
      .flatten()
      .buffer(size:2,skip:1)
      .flatten()
      .mix(full_coverage_bigwigs)
      .collect()
      .flatten()
      .toSortedList()
      .set{ full_coverage_inputs } */
  
    /*
     * Step 7b: Create Full Coverage Objects
     */

/*     process sampleCreateCoverageObjects {

        echo true
        tag "Creating Coverage Objects [ $full_coverage_input ]"
        publishDir "${params.basedir}/Coverage_Objects",mode:'copy'

        input:
        file fullCov_file from fullCov_file
        file fullCov_samples_manifest from fullCov_samples_manifest
        file full_coverage_input from full_coverage_inputs
        file inferred_strand_R_object from inferred_strand_objects

        output:
        file "*"

        shell:
        if (params.sample == "paired") {
            coverage_pe = "TRUE"
        }
        if (params.sample == "single") {
            coverage_pe = "FALSE"
        }
        coverage_cores = "${params.coverage_cores}"
        coverage_reference = "${params.reference}"
        coverage_experiment = "${params.experiment}"
        coverage_prefix = "${params.experiment_prefix}"
        coverage_fullCov = "TRUE"
        '''
        Rscript !{fullCov_file} -o !{coverage_reference} -m . -e !{coverage_experiment} -p !{coverage_prefix} -l !{coverage_pe} -f !{coverage_fullCov} -c !{coverage_cores}
        '''
    }
} */

/* if (params.step8) { */

    /*
     * Step 8: Call Variants
     */

/*     variant_calls_bam
      .combine(snvbed)
      .combine(variant_assembly)
      .set{ variant_calls }

    process sampleVariantCalls {

        echo true
        tag "Prefix: $variant_bams_prefix | Sample: [ $variant_calls_bam_file, $variant_calls_bai ]"
        publishDir "${params.basedir}/Variant_Calls",mode:'copy'

        input:
        set val(variant_bams_prefix), file(variant_calls_bam_file), file(variant_calls_bai), file(snv_bed), file(variant_assembly_file) from variant_calls

        output:
        file "${variant_bams_prefix}.vcf.gz" into compressed_variant_calls
        file "${variant_bams_prefix}.vcf.gz.tbi" into compressed_variant_calls_tbi

        shell:
        snptmp = "${variant_bams_prefix}_tmp.vcf"
        snpoutgz = "${variant_bams_prefix}.vcf.gz"
        '''
        samtools mpileup \
        -l !{snv_bed} \
        -AB \
        -q0 \
        -Q13 \
        -d1000000 \
        -uf !{variant_assembly_file} !{variant_calls_bam_file} \
        -o !{snptmp}
        bcftools call \
        -mv \
        -Oz !{snptmp} > !{snpoutgz}
        tabix -p vcf !{snpoutgz}
        '''
    } */

/*     compressed_variant_calls
      .flatten()
      .collect()
      .toSortedList()
      .set{ collected_variant_calls }

    compressed_variant_calls_tbi
      .flatten()
      .collect()
      .toSortedList()
      .set{ collected_variant_calls_tbi } */


    /*
     * Step 8b: Merge Variant Calls
     */

/*     process sampleVariantCallsMerge {

        echo true
        tag "Samples: $collected_variants"
        publishDir "${params.basedir}/Merged_Variants",mode:'copy'

        input:
        file collected_variants from collected_variant_calls
        file collected_variants_tbi from collected_variant_calls_tbi

        output:
        file "*"

        shell:
        '''
        vcf-merge !{collected_variants} | bgzip -c > mergedVariants.vcf.gz
        '''
    }
} */

/*
 * Step 9: Expressed Regions
 */

/* process sampleExpressedRegions {

    echo true
    tag "Sample: $expressed_regions_mean_bigwig"
    publishDir "${params.basedir}/Expressed_Regions",mode:'copy'

    input:
    file expressedRegions_file from expressedRegions_file
    file chr_sizes from chr_sizes
    file expressed_regions_mean_bigwig from expressed_regions_mean_bigwigs

    output:
    file "*" 

    shell:
    expressed_regions_cores = "${params.expressedregion_cores}"
    '''
    Rscript !{expressedRegions_file} \
    -m !{expressed_regions_mean_bigwig} \
    -o . \
    -i !{chr_sizes} \
    -c !{expressed_regions_cores}
    '''
}
 */
