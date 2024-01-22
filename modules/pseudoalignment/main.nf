// This module contains processes for pseudoalignment to a reference
// transcriptome (including processes for both Salmon and Kallisto)

process TxQuantSalmon {

    tag "$prefix"
    publishDir "${params.output}/salmon_tx/${prefix}",mode:'copy'

    input:
        path salmon_index
        path complete_manifest
        tuple val(prefix), path(salmon_inputs)

    output:
        path "${prefix}/*"
        path "${prefix}_quant.sf", emit: tx_quants
        path "tx_quant_${prefix}.log"

    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness and determine strand flag
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            strand_flag="SF"
        elif [ $strand == 'reverse' ]; then
            strand_flag="SR"
        else
            strand_flag="U"
        fi
        
        if [ !{params.sample} == "paired" ]; then
            strand_flag="I$strand_flag"
            sample_flag="-1 !{prefix}*trimmed*_1.f*q* -2 !{prefix}*trimmed*_2.f*q*"
        else
            sample_flag="-r !{prefix}*.f*q*"
        fi
        
        !{params.salmon} quant \
            -i $PWD \
            -p !{task.cpus} \
            -l ${strand_flag} \
            ${sample_flag} \
            -o !{prefix} \
            !{params.salmon_quant_args}

        cp !{prefix}/quant.sf !{prefix}_quant.sf
        cp .command.log tx_quant_!{prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

process TxQuantKallisto {

    tag "$prefix"
    publishDir "${params.output}/kallisto_tx/${prefix}", mode:'copy'

    input:
        path kallisto_index
        path complete_manifest
        tuple val(prefix), path(fastqs)

    output:
        path "abundance.h5"
        path "run_info.json"
        path "tx_quant_${prefix}.log"
        path "${prefix}_abundance.tsv", emit: tx_quants
        
    shell:
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Find this sample's strandness
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            strand_flag="--fr-stranded"
        elif [ $strand == 'reverse' ]; then
            strand_flag="--rf-stranded"
        else
            strand_flag=""
        fi
        
        #  Run the quantification step
        if [ !{params.sample} == "paired" ]; then
            !{params.kallisto} quant \
                -t !{task.cpus} \
                $strand_flag \
                -i !{kallisto_index} \
                -o . \
                !{params.kallisto_quant_paired_args} \
                !{prefix}*trimmed*_1.f*q* !{prefix}*trimmed*_2.f*q*
        else
            !{params.kallisto} quant \
                -t !{task.cpus} \
                !{params.kallisto_quant_single_args} \
                $strand_flag \
                -i !{kallisto_index} \
                -o . \
                *.f*q*
        fi
        
        mv abundance.tsv !{prefix}_abundance.tsv
        cp .command.log tx_quant_!{prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}
