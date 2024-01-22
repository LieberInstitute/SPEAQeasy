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
