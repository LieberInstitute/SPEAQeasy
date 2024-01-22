// This module includes processes for calling variants at a small list of highly
// variable SNVs, a step performed for human data, and intended to allow
// downstream resolution of potential sample-id issues (see
// http://research.libd.org/SPEAQeasy-example/swap_speaqeasy.html)

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


// Merge Variant Calls
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
