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
