// This module includes processes related to quantifying coverage and finding
// expressed regions. These processes may be enabled with the '--coverage' and
// '--fullCov' options.

process Coverage {

    tag "$prefix"
    publishDir "${params.output}/coverage/wigs", mode:'copy'

    input:
        path complete_manifest
        tuple val(prefix), path(sorted_bam), path(sorted_bam_index)
        path chr_sizes

    output:
        path "${prefix}*.wig", emit: wig_files
        path "bam2wig_${prefix}.log"

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            if [ !{params.sample} == "paired" ]; then
                strand_flag='-d 1++,1--,2+-,2-+'
            else
                strand_flag='-d ++,--'
            fi
        elif [ $strand == 'reverse' ]; then
            if [ !{params.sample} == "paired" ]; then
                strand_flag='-d 1+-,1-+,2++,2--'
            else
                strand_flag='-d +-,-+'
            fi
        else
            strand_flag=""
        fi

        python3 $(which bam2wig.py) \
            -s !{chr_sizes} \
            -i !{sorted_bam} \
            !{params.bam2wig_args} \
            -o !{prefix} \
            $strand_flag

        cp .command.log bam2wig_!{prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2-|| true)
        echo "$temp" > bash_vars.txt
        '''
}

// Convert wig files to bigWig format    
process WigToBigWig {

    tag "$prefix"
    publishDir "${params.output}/coverage/bigWigs", mode: 'copy'
    
    input:
        tuple val(prefix), path(wig_file)
        path chr_sizes
    
    output:
        path "*.bw", emit: coverage_bigwigs
    
    shell:
        '''
        !{params.wigToBigWig} \
            !{params.wigToBigWig_args} \
            !{wig_file} \
            !{chr_sizes} \
            !{prefix}.bw
        '''
}

// Compute mean coverage across samples by strand
process MeanCoverage {

    tag "Strand: ${read_type}"
    publishDir "${params.output}/coverage/mean", mode: 'copy'

    input:
        tuple val(read_type), path(mean_coverage_bigwig)
        path chr_sizes
        path chunk_apply_script

    output:
        tuple val(read_type), path("mean*.bw"), emit: mean_bigwigs
        path "mean_coverage_${read_type}.log"

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        if [ !{read_type} == "unstranded" ] ; then
            outwig="mean"
        elif [ !{read_type} == "forward" ]; then
            outwig="mean.forward"
        else
            outwig="mean.reverse"
        fi
        
        precision=8
        
        file_regex='.*\\.(bw|wig)$'
        command="!{params.wiggletools} write temp_wig[index].wig sum [files]"

        #   To compute a mean, we'll sum all files then multiply by the
        #   scale factor (reciprocal of number of files)
        num_files=$(ls -1 | grep -E "$file_regex" | wc -l)
        scale_factor=$(!{params.bc} <<< "scale=${precision};1/$num_files")
        echo "Scale factor is $scale_factor."
        
        #   Break the summation command into chunks (necessary for large
        #   datasets where the number of open file handles may be exceeded)
        bash !{chunk_apply_script} !{params.wiggletools_batch_size} "$file_regex" "$command"

        #   Note that even if one file is left, this command is necessary
        #   because of the scaling (mean instead of sum)
        echo "Computing a mean from the final files..."
        file_list=$(ls -1 | grep -E "$file_regex")
        !{params.wiggletools} write $outwig.wig scale $scale_factor sum $file_list
        echo "Done."
        
        !{params.wigToBigWig} $outwig.wig !{chr_sizes} $outwig.bw
        
        cp .command.log mean_coverage_!{read_type}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

process CoverageObjects {
    
    tag "Strand: ${read_type}"
    publishDir "${params.output}/coverage_objects", mode: 'copy'

    input:
        path fullCov_script
        path complete_manifest
        tuple val(read_type), path(bigwigs)

    output:
        file "*"

    shell:
        '''
        !{params.Rscript} !{fullCov_script} \
            -o !{params.reference} \
            -e !{params.experiment} \
            -c !{task.cpus} \
            -s !{read_type}
        
        cp .command.log coverage_objects.log
        '''
}

process ExpressedRegions {
    
    tag "Strand: ${read_type}"
    publishDir "${params.output}/expressed_regions", mode:'copy'
    
    input:
        path expressed_regions_script
        path chr_sizes
        tuple val(read_type), path(mean_bigwigs)
    
    output:
        file "*"
    
    shell:
        '''
        for meanfile in ./mean*.bw; do
            !{params.Rscript} !{expressed_regions_script} \
                -m ${meanfile} \
                -o . \
                -i !{chr_sizes} \
                -c !{task.cpus}
        done
        
        cp .command.log expressed_regions_!{read_type}.log
        '''
}
