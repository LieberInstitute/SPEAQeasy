// This module contains all processes performing alignment of raw FASTQ reads to
// the reference genome, involving either HISAT2 or STAR depending on settings

/*
 * Alignment (HISAT2 by default, or otherwise STAR)
 */

process AlignStar {
    tag "$prefix"
    publishDir "${params.output}/alignment", mode:'copy'
    
    input:
        path star_index
        tuple val(prefix), path(single_star_input)
        
    output:
        tuple val(prefix), path("*.bam"), emit: alignment_output
        path "${prefix}_unmapped_*.fastq", optional: true
        path "*_STAR_alignment.log", emit: alignment_summaries
        
    shell:            
        if (params.sample == "paired" && params.unalign) {
            unaligned_args = "--outReadsUnmapped Fastx"
        } else {
            unaligned_args = ""
        }
        '''
        #  Recreate the index directory (nextflow cannot handle directories
        #  as channel input). This step is likely sensitive to the STAR
        #  version used
        mkdir index_dir
        mv *.txt index_dir/
        mv *.tab index_dir/
        mv SA index_dir/
        mv SAindex index_dir/
        mv Genome index_dir/
        
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Determine file names for input FASTQs
        if [ !{params.sample} == "paired" ]; then
            fastq_files="!{prefix}*trimmed*_1.f*q* !{prefix}*trimmed*_2.f*q*"
        else
            fastq_files="*.f*q*"
        fi
        
        #  Tell STAR to decompress files if needed 
        if [ "$(echo $fastq_files | rev | cut -d "." -f 1 | rev)" == "gz" ]; then
            compress_args='--readFilesCommand gunzip -c'
        else
            compress_args=''
        fi
        
        #  Perform alignment
        !{params.star} \
            --genomeDir ./index_dir \
            --runThreadN !{task.cpus} \
            --readFilesIn ${fastq_files} \
            --outSAMtype BAM Unsorted \
            ${compress_args} \
            !{unaligned_args} \
            !{params.star_args}
            
        #  Adjust file names (make unique)
        mv Log.final.out !{prefix}_STAR_alignment.log
        mv Aligned.out.bam !{prefix}.bam
        if [ -f Unmapped.out.mate1 ]; then
            mv Unmapped.out.mate1 !{prefix}_unmapped_mate1.fastq
            mv Unmapped.out.mate2 !{prefix}_unmapped_mate2.fastq
        fi
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}
        
process SingleEndHisat {

    tag "$prefix"
    publishDir "${params.output}/alignment",mode:'copy'

    input:
        path hisat_index
        path complete_manifest
        tuple val(prefix), path(single_hisat_input)

    output:
        tuple val(prefix), path("${prefix}.bam"), emit: alignment_output
        path "*_align_summary.txt", emit: alignment_summaries

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ ${strand} == "unstranded" ]; then
            hisat_strand=""
        elif [ ${strand} == "forward" ]; then
            hisat_strand="--rna-strandness F"
        else
            hisat_strand="--rna-strandness R"
        fi
        
        #  Run Hisat2
        !{params.hisat2} \
            -p !{task.cpus} \
            -x !{params.annotation}/reference/!{params.reference}/assembly/index/hisat2_assembly_!{params.anno_suffix} \
            -U !{single_hisat_input} \
            ${hisat_strand} \
            !{params.hisat2_args} \
            2> !{prefix}_align_summary.txt \
        | !{params.samtools} view -b -F 4 -o !{prefix}.bam
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

process PairedEndHisat {

    tag "$prefix"
    publishDir "${params.output}/alignment",'mode':'copy'

    input:
        path hisat_index
        path complete_manifest
        tuple val(prefix), path(fq_files)

    output:
        tuple val(prefix), path("${prefix}.bam"), emit: alignment_output
        path "*_unpaired*.fastq" optional true
        path "*_align_summary.txt", emit: alignment_summaries

    shell:
        if (params.unalign) {
            unaligned_opt = "--un-conc ${prefix}_discordant.fastq"
        } else {
            unaligned_opt = ""
        }
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ ${strand} == "unstranded" ]; then
            hisat_strand=""
        elif [ ${strand} == "forward" ]; then
            hisat_strand="--rna-strandness FR"
        else
            hisat_strand="--rna-strandness RF"
        fi
        
        #  If this sample had unpaired trimming outputs, include them
        if [ "!{params.keep_unpaired}" == "true" ]; then
            unpaired_opt='-U !{prefix}_unpaired_1.fastq,!{prefix}_unpaired_2.fastq'
        else
            unpaired_opt=''
        fi
        
        #  Run Hisat2
        !{params.hisat2} \
            -p !{task.cpus} \
            -x !{params.annotation}/reference/!{params.reference}/assembly/index/hisat2_assembly_!{params.anno_suffix} \
            -1 !{prefix}*trimmed*_1.f*q* \
            -2 !{prefix}*trimmed*_2.f*q* \
            ${unpaired_opt} \
            ${hisat_strand} \
            !{params.hisat2_args} \
            !{unaligned_opt} \
            2> !{prefix}_align_summary.txt \
        | !{params.samtools} view -b -F 4 -o !{prefix}.bam
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

/*
 * Sort and index the aligned BAM
 */

process BamSort {
	
    tag "$prefix"
    publishDir "${params.output}/alignment/bam_sort",'mode':'copy'

    input:
        tuple val(prefix), path(input_bam)

    output:
        tuple val(prefix), path("${prefix}_sorted.bam"), path("${prefix}*_sorted.bam.bai"), emit: sorted_bams

    shell:
        '''
        !{params.samtools} sort -T temporary -l 6 --no-PG -@ !{task.cpus} !{input_bam} -o !{prefix}_sorted.bam
        !{params.samtools} index !{prefix}_sorted.bam
        '''
}
