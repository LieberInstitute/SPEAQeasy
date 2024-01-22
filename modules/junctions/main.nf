// This module contains processes involved in quantifying junctions after
// allignment to the reference genome

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
