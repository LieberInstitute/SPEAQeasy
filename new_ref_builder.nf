#!/usr/bin/env nextflow

/*
Pipeline Overview:

        -   _000a_hisat2_index_builder

        -  IIa: Download GENCODE GTF File
        -  IIb: Build Bed File
        - IIIa: Download Salmon TX FA
        - IIIb: Build Salmon Index
*/

/* WinterFlow new variables */

// Define Configurable Variables
//params.annotations = false
params.annotations = "./Annotation"
params.reference = false
params.output = false

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

// Annotation Path Validation
// disabled since does not seem to be asigning values by default...
// if (!params.annotations) {
//         params.annotations = "./Annotation"
// }

params.cores = '4'
params.ercc_cores = '4'
params.salmon_cores = '4'
params.hisat_cores = '4'

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define External Scripts
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//##TODO(iaguilar): This should live in an external config file... (Dev ######)
// prep_bed = file("${params.scripts}/prep_bed.R")
// check_R_packages_script = file("${params.scripts}/check_R_packages.R")

if (!params.output) {
	params.basedir = "${params.output}"
	params.index_out = "${params.annotations}"
}

if (params.reference == "hg38") {

	// Step 3: hisat2
	params.fa_link = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz"
/* commented since the WineerFlow update should resolve redundancies */
	// params.fa_gz = "GRCh38.primary_assembly.genome.fa.gz"
	// params.fa = "GRCh38.primary_assembly.genome.fa"
	// params.hisat_prefix = "hisat2_GRCh38primary"
	// params.hisat_assembly = "GENCODE/GRCh38_hg38/assembly"

	// Step 4: gencode gtf
//##TODO(iaguilar): Check if gtf_link and gtf_gz are not redundant since link includes fa_gaz value (Dev ######)
	params.gencode_gtf_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"
	params.gencode_gtf_gz = "gencode.v25.annotation.gtf.gz"
	params.gencode_gtf = "gencode.v25.annotation.gtf"
	params.feature_output_prefix = "Gencode.v25.hg38"

	// Step 5: python coverage
//##TODO(iaguilar): Briefly explain what is this file used for (Doc ######)
	//chr_sizes = file("${params.annotations}/chrom_sizes/hg38.chrom.sizes.gencode")

	// Step 6: salmon
//##TODO(iaguilar): Explain why step 6 is enabled if reference is hg38...  (Doc ######)
	params.step6 = true
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
	params.tx_fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.transcripts.fa.gz"
	params.tx_fa_gz = "gencode.v25.transcripts.fa.gz"
	params.tx_fa = "gencode.v25.transcripts.fa"
	params.salmon_prefix = "salmon_0.8.2_index_gencode.v25.transcripts"
	params.salmon_assembly = "GENCODE/GRCh38_hg38/transcripts"
}

if (params.reference == "hg19") {

	// Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
	params.fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
/* commented since the WineerFlow update should resolve redundancies */
  // params.fa_gz = "GRCh37.primary_assembly.genome.fa.gz"
	// params.fa = "GRCh37.primary_assembly.genome.fa"
	// params.hisat_prefix = "hisat2_GRCh37primary"
	// params.hisat_assembly = "GENCODE/GRCh37_hg19/assembly"

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

}
if (params.reference == "mm10") {

	// Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gaz value (Dev ######)
	params.fa_link = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/GRCm38.primary_assembly.genome.fa.gz"
/* commented since the WineerFlow update should resolve redundancies */
  // params.fa_gz = "GRCm38.primary_assembly.genome.fa.gz"
	// params.fa = "GRCm38.primary_assembly.genome.fa"
	// params.hisat_prefix = "GRCm38_mmhisat2_GRCm38primary"
	// params.hisat_assembly = "GENCODE/GRCm38_mm10/assembly"

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

}
if (params.reference == "rn6") {

	// Step 3: hisat2
//##TODO(iaguilar): Check if fa_link and fa_gz are not redundant since link includes fa_gz value (Dev ######)
	params.fa_link = "ftp://ftp.ensembl.org/pub/release-86/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
/* commented since the WineerFlow update should resolve redundancies */
  // params.fa_gz = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz"
	// params.fa = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
	// params.hisat_prefix = "hisat2_Rnor6.0toplevel"
	// params.hisat_assembly = "ensembl/Rnor_6.0"

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

}

/*
  * DEFINE PATHS TO MK MODULES
*/

/* _000a_hisat2_index_builder */
module_mk_000a_hisat2_index_builder = "./mkmodules/reference_building/hisat2_index_builder//"

/*
  * BEGIN PIPELINE
*/

/* Read mkfile module files */
Channel
        .fromPath("${module_mk_000a_hisat2_index_builder}/*")
        .toList()
        .set{ mkfiles_000a }

process _000a_hisat2_index_builder {

	tag "Building ht2 index at: ${params.annotations}/hisat2/${params.reference}/"

  /* The storeDir directive will be used to define the output directory as a cache for hisat2 references
   * Name of the cache dir will be dynamically constructed
  */
  storeDir "${params.annotations}/hisat2/${params.reference}/"

  input:
  file mk_file from mkfiles_000a

	output:
	file "*.ht2" into ht2_index_from_000a_hisat2_index_builder
  file "*.fa*" into fasta_from_000a_hisat2_index_builder

	"""
  bash runmk.sh REFERENCE_FASTA_URL="${params.fa_link}" NUMBER_OF_THREADS="${params.hisat_cores}"
	"""
}
