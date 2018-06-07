#!/usr/bin/env nextflow

/*
Pipeline Overview:

   Preprocessing:
        -   Ia: Download Assembly FA
        -   Ib: Build HISAT Index
        -  IIa: Download GENCODE GTF File
        -  IIb: Build Bed File
        - IIIa: Download Salmon TX FA
        - IIIb: Build Salmon Index
*/

// Define Configurable Variables
params.annotation = false
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
if (!params.annotation) {
        params.annotations = "./Annotation"
}

params.cores = '4'
params.ercc_cores = '4'
params.salmon_cores = '4'
params.hisat_cores = '4'

// External Script Path Validation
//##TODO(iaguilar): This param was not defined neither in help nor in the variable definition block (Dev ######)
//##TODO(iaguilar): Comment in original dev was: "It's the directory where the shell files are located at. You only need to specify it if you cloned this repository somewhere else"; since this scripts folder will be versioned with the NF pipeline, there is no need to allow it to be on another directory (Dev ######)
params.scripts = "./scripts"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define External Scripts
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//##TODO(iaguilar): This should live in an external config file... (Dev ######)
prep_bed = file("${params.scripts}/prep_bed.R")
check_R_packages_script = file("${params.scripts}/check_R_packages.R")

if (!params.output) {
	params.basedir = "${params.output}"
	params.index_out = "${params.annotations}"
}

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

}


// BEGIN PIPELINE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Step Ia: GENCODE Assembly FA
 */

/* ###################################
	STARTS REFERENCES BUILDING BLOCK 
################################### */

/*
 * Step Ia: GENCODE Assembly FA
 */

// Define the output directory path for Genome Reference dir, and hisat2 index dir, based on previously defined parameters
params.hisat_idx_output = "${params.index_out}/${params.hisat_assembly}"

if ( ! file("${params.hisat_idx_output}/fa/${params.fa}").exists() ){ 
//If file does not exists, launch process to create it
	process pullGENCODEassemblyfa {
		
	tag "Downloading Assembly FA File: ${params.fa}"
	publishDir "${params.hisat_idx_output}/fa",'mode':'copy'

	output:
	file("${params.fa}") into reference_assembly
	file("${params.fa}") into variant_assembly

	script:
	fa_gz_file_link = "${params.fa_link}"
	fa_gz = "${params.fa_gz}"
	"""
		wget "$fa_gz_file_link"
		gunzip "$fa_gz"
	"""
		}
} else {
		//If file already exists, load it from path into the corresponding channel
		println "[WG-LOG] Skipping download process for ${params.hisat_idx_output}/fa/${params.fa}, file already exist"
		Channel
			.fromPath("${params.hisat_idx_output}/fa/${params.fa}")
			.into{ reference_assembly; variant_assembly }
}

/*
 * Step Ib: Build HISAT Index
 */

/*Check if the Hisat2 index directory exists*/
if ( ! file("${params.hisat_idx_output}/index").exists() ) {
/* If the index directory does not exist, launch process to create them */
	println "[WG-LOG] building hisat2 index at ${params.hisat_idx_output}/index/"
	process buildHISATindex {

		
		tag "Building HISAT2 Index: ${params.hisat_prefix}"
		publishDir "${params.hisat_idx_output}/index",'mode':'copy'

		input:
		file reference_fasta from reference_assembly

		output:
		file("${params.hisat_prefix}.*") into hisat_index_built

		script:
		prefix = "${params.hisat_prefix}"
		/* $prefix is "the hisat2_index_base. Write ht2 data to files with this dir/basename" */
		"""
			${params.hisat2build} -p "${params.hisat_cores}" $reference_fasta $prefix
		"""
	}
} else {
/* If index already exists, load files from path into the corresponding channel */
	println "[WG-LOG] Skipping hisat2 index process for ${params.hisat_idx_output}/index/, index files already exist"
	Channel
		.fromPath("${params.hisat_idx_output}/index/${params.hisat_prefix}.*.ht2")
		.set{ hisat_index_built }
}

/* Channel post-processing */
hisat_index_built // get every *.ht2 file in this channel
	.toSortedList() // group every *.ht2 item into a single, sorted list
	// /* .dump() */ //For debug, make channel content informative if you use the -dump-channels option of NF
	.set{ hisat_index } // pass *.ht2 list to a new channel

/*
 * Step IIa: GENCODE GTF Download
 */

// Define the output directory path for GTF file, and the prep bed, based on previously defined parameters
params.gencode_gtf_out = "${params.index_out}/RSeQC/${params.reference}"

if ( ! file("${params.gencode_gtf_out}/gtf/${params.gencode_gtf}").exists() ){
//If file does not exists, launch process to create it
println "[WG-LOG] downloading ${params.gencode_gtf_out}/gtf/${params.gencode_gtf}"
	process pullGENCODEgtf {

		
		tag "Downloading GTF File: ${params.gencode_gtf}"
		publishDir "${params.gencode_gtf_out}/gtf",'mode':'copy'

		output:
		file("${params.gencode_gtf}") into gencode_gtf
		file("${params.gencode_gtf}") into create_counts_gtf
		file("${params.gencode_gtf}") into gencode_feature_gtf

		script:
		gencode_gtf_link = "${params.gencode_gtf_link}"
		gencode_gtf_gz = "${params.gencode_gtf_gz}"
		gencode_gtf_file = "${params.gencode_gtf}"
		"""
			wget "$gencode_gtf_link"
			gunzip "$gencode_gtf_gz"
		"""
	}
} else {
	//If file already exists, load it from path into the corresponding channel
	println "[WG-LOG] Skipping download process for ${params.gencode_gtf_out}/gtf/${params.gencode_gtf}, file already exists"
	Channel
		.fromPath("${params.gencode_gtf_out}/gtf/${params.gencode_gtf}")
		.into{ gencode_gtf; create_counts_gtf; gencode_feature_gtf }
}


/*
 * Step IIb: Build Bed File
 */

if ( ! file("${params.gencode_gtf_out}/bed/${params.reference}.bed").exists() ){
//If file does not exists, launch process to create it
println "[WG-LOG] building ${params.gencode_gtf_out}/bed/${params.reference}.bed"
	process buildPrepBED {

		
		tag "Building Bed File: ${params.reference}"
		publishDir "${params.gencode_gtf_out}/bed",'mode':'copy'

		input:
		file gencode_gtf from gencode_gtf
		file prep_bed from prep_bed
		file check_R_packages_script from check_R_packages_script

		output:
		file("${name}.bed") into bedfile

		shell:
		name = "${params.reference}"
		'''
		Rscript !{check_R_packages_script} \
		&& Rscript !{prep_bed} -f !{gencode_gtf} -n !{name}
		'''
	}
} else {
		//If file already exists, load it from path into the corresponding channel
		println "[WG-LOG] Skipping build process for ${params.gencode_gtf_out}/bed/${params.reference}.bed, file already exists"
		Channel
			.fromPath("${params.gencode_gtf_out}/bed/${params.reference}.bed")
			.set{ bedfile }
}



// for hg38, hg19, and mm10, step 6 is enabled by params.step6 = true 
// during Define Reference Paths/Scripts + Reference Dependent Parameters
if (params.step6) {

	params.salmon_idx_output = "${params.index_out}/${params.salmon_assembly}"

	/*
	 * Step IIIa: GENCODE TX FA Download
	 */
	// use the file operator to find if a file exists in the system
	if ( ! file("${params.salmon_idx_output}/fa/${params.tx_fa}").exists() ) {
		//If file does not exists, launch process for download
		println "[WG-LOG] downloading ${params.salmon_idx_output}/fa/${params.tx_fa}"
		process pullGENCODEtranscripts {

			
			tag "Downloading TX FA File: ${params.tx_fa}"
			publishDir "${params.salmon_idx_output}/fa",'mode':'copy'

			output:
			file("${params.tx_fa}") into transcript_fa
			//TODO (iaguilar) put shell code in same style as for fasta reference download (Dev ###)
			script:
			tx_fa_link = "${params.tx_fa_link}"
			tx_fa_gz = "${params.tx_fa_gz}"
			"""
				wget $tx_fa_link
				gunzip $tx_fa_gz
			"""
		}
	} else {
		//If file already exists, load it from path into the corresponding channel
		println "[WG-LOG] Skipping download process for ${params.salmon_idx_output}/fa/${params.tx_fa}, file already exists"
		Channel
			.fromPath("${params.salmon_idx_output}/fa/${params.tx_fa}")
			.set{ transcript_fa }
	}

	/*
	 * Step IIIb: Salmon Transcript Build
	 */

/* NOT FULLY WORKING; seems the salmon_index channel manipulation in this version affects the TXQUANT process; only one sample gets processed
if ( ! file("${params.salmon_idx_output}/salmon/${params.salmon_prefix}").exists() ){
		//If file does not exists, launch process to create it
		println "[WG-LOG] creating ${params.salmon_idx_output}/salmon/${params.salmon_prefix}/*"
		process buildSALMONindex {
			
			tag "Building Salmon Index: ${params.salmon_prefix}"
			publishDir "${params.salmon_idx_output}/salmon",'mode':'copy'
			input:
			file tx_file from transcript_fa
			output:
			file("${params.salmon_prefix}") into salmon_index
			script:
			salmon_idx = "${params.salmon_prefix}"
			"""
				${params.salmon} index -t $tx_file -i $salmon_idx -p ${params.salmon_cores} --type quasi -k 31
			"""
		}
	} else {
		//If directory already exists, load it from path into the corresponding channel
		println "[WG-LOG] Skipping build process for ${params.salmon_idx_output}/salmon/${params.salmon_prefix}, directory already exist"
		Channel
			.fromPath("${params.salmon_idx_output}/salmon/${params.salmon_prefix}")
			.set{ salmon_index }
	}
} */

/* Using the old way of triggering reference constructions */
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
		"""
			${params.salmon} index -t $tx_file -i $salmon_idx -p ${params.salmon_cores} --type quasi -k 31
		"""
	}

	// Read the salmon index from path, and/or mix with the channel output from the buildSALMONindex process
	// This effectively means that whether the index was created or it was already present, the flow can continue
	 Channel
		.fromPath("${params.salmon_idx_output}/salmon/${params.salmon_prefix}")
		.mix(salmon_index_built)
		.toSortedList()
		.flatten()
		.distinct()
		.toSortedList()
		.set{ salmon_index }
}
