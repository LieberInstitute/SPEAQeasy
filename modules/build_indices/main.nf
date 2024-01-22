// This module contains all processes related to building an index for a
// software tool. STAR, HISAT2, Salmon, and Kallisto all need indices
// built for the particular reference files/ settings chosen; the associated
// processes are defined here. 'BuildAnnotationObjects' is also defined here,
// which is a custom R script to generate objects used in the 'CountObjects'
// process

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
