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
