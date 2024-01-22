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
