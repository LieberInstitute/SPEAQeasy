// This module contains QC-related processes that occur early in the SPEAQeasy
// workflow. In particular, it includes 'Trimming', FastQC before and after
// trimming, and ERCC quantification

def get_file_ext(f) {
    if (f.name.toString().tokenize(".")[-1] == "gz") {
        return('.fastq.gz')
    } else {
        return('.fastq')
    }
}

/*
 * Perform FastQC on the untrimmed reads, as an initial quality gauge
 */

process QualityUntrimmed {

    tag "$prefix"
    publishDir "${params.output}/fastQC/untrimmed", mode:'copy', pattern:'*_fastqc'

    input:
        tuple val(prefix), path(fastqc_untrimmed_input)

    output:
        path "${prefix}*_fastqc"
        tuple val(prefix), path("*_summary.txt"), emit: quality_reports
        path "*_fastqc_data.txt", emit: count_objects_quality_metrics_untrimmed

    shell:
        if (params.sample == "single") {
            copy_command = "cp ${prefix}_fastqc/summary.txt ${prefix}_untrimmed_summary.txt"
            data_command = "cp ${prefix}_fastqc/fastqc_data.txt ${prefix}_untrimmed_fastqc_data.txt"
        } else {
            copy_command = "cp ${prefix}_1_fastqc/summary.txt ${prefix}_1_untrimmed_summary.txt && cp ${prefix}_2_fastqc/summary.txt ${prefix}_2_untrimmed_summary.txt"
            data_command = "cp ${prefix}_1_fastqc/fastqc_data.txt ${prefix}_1_untrimmed_fastqc_data.txt && cp ${prefix}_2_fastqc/fastqc_data.txt ${prefix}_2_untrimmed_fastqc_data.txt"
        }
        '''
        !{params.fastqc} -t !{task.cpus} --extract !{params.fastqc_args} !{fastqc_untrimmed_input}
        !{copy_command}
        !{data_command}
        '''
}

/*
 * Trim samples (dependent on user-chosen settings)
 */
 
process Trimming {

    tag "$fq_prefix"
    publishDir "${params.output}/trimming", mode:'copy', pattern:'*_trimmed*.f*q{.gz,}'

    input:
        tuple val(fq_prefix), path(fq_summary), path(fq_file)

    output:
        tuple val(fq_prefix), path("${fq_prefix}_trimmed*.f*q{.gz,}"), optional: true, emit: trimmed_fastqc_inputs_raw
        tuple val(fq_prefix), path("${fq_prefix}*.f*q{.gz,}"), emit: trimming_outputs

    shell:
        file_ext = get_file_ext(fq_file[0])
        if (params.sample == "single") {
            output_option = "${fq_prefix}_trimmed.fastq"
            trim_mode = "SE"
            adapter_fa_temp = params.adapter_fasta_single
            trim_clip = params.trim_adapter_args_single
        } else {
            output_option = "${fq_prefix}_trimmed_paired_1.fastq ${fq_prefix}_unpaired_1.fastq ${fq_prefix}_trimmed_paired_2.fastq ${fq_prefix}_unpaired_2.fastq"
            trim_mode = "PE"
            adapter_fa_temp = params.adapter_fasta_paired
            trim_clip = params.trim_adapter_args_paired
        }
        '''
        ( set -o posix ; set ) > bash_vars.txt
        
        #  Determine whether to trim the FASTQ file(s). This is ultimately
        #  controlled by the '--trim_mode' command flag.
        if [ "!{params.trim_mode}" == "force" ]; then
            do_trim=true
        elif [ "!{params.trim_mode}" == "skip" ]; then
            do_trim=false
        elif [ "!{params.sample}" == "single" ]; then
            #  Then '--trim_mode "adaptive"' was selected, and data is single-end
            if [ $(grep "Adapter Content" !{fq_summary} | cut -f 1)  == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        else
            #  Then '--trim_mode "adaptive"' was selected, and data is paired-end
            result1=$(grep "Adapter Content" !{fq_prefix}_1_untrimmed_summary.txt | cut -c1-4)
            result2=$(grep "Adapter Content" !{fq_prefix}_2_untrimmed_summary.txt | cut -c1-4)
            if [ $result1 == "FAIL" ] || [ $result2 == "FAIL" ]; then
                do_trim=true
            else
                do_trim=false
            fi
        fi
        
        #  Run trimming if required
        if [ "$do_trim" == true ]; then
            #  This solves the problem of trimmomatic and the adapter fasta
            #  needing hard paths, even when on the PATH.
            if [ !{params.use_long_paths} == "true" ]; then
                trim_jar=!{params.trimmomatic}
                adapter_fa=!{adapter_fa_temp}
            else
                trim_jar=$(which !{params.trimmomatic})
                adapter_fa=$(which !{adapter_fa_temp})
            fi
            
            #  Now determine the arguments to pass to trimmomatic regarding
            #  adapter trimming. The flexibility here allows the user to 
            #  potentially not perform adapter trimming.
            if [ "!{trim_clip}" == "" ]; then
                adapter_trim_settings=""
            else
                adapter_trim_settings="ILLUMINACLIP:$adapter_fa:!{trim_clip}"
            fi

            java -Xmx512M \
                -jar $trim_jar \
                !{trim_mode} \
                -threads !{task.cpus} \
                -phred33 \
                *.f*q* \
                !{output_option} \
                $adapter_trim_settings \
                !{params.trim_quality_args}
        else
            #  Otherwise rename files (signal to nextflow to output these files)
            if [ "!{params.sample}" == "single" ]; then
                mv !{fq_prefix}!{file_ext} !{fq_prefix}_untrimmed!{file_ext}
            else
                mv !{fq_prefix}_1!{file_ext} !{fq_prefix}_untrimmed_1!{file_ext}
                mv !{fq_prefix}_2!{file_ext} !{fq_prefix}_untrimmed_2!{file_ext}
            fi
        fi
        
        cp .command.log trimming_!{fq_prefix}.log
        
        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep ">" | cut -d " " -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}

/*
 * Run FastQC Quality Check on Trimmed Files
 */

process QualityTrimmed {

    tag "$prefix"
    publishDir "${params.output}/fastQC/trimmed",'mode':'copy', pattern:'*_fastqc'

    input:
        tuple val(prefix), path(fastqc_trimmed_input)

    output:
        path "${prefix}*_fastqc"
        path "${prefix}*_trimmed_*.txt", emit: trimmed_fastqc_outs

    shell:
        if (params.sample == "single") {
            copy_command = "cp ${prefix}_trimmed_fastqc/summary.txt ${prefix}_trimmed_summary.txt"
            data_command = "cp ${prefix}_trimmed_fastqc/fastqc_data.txt ${prefix}_trimmed_fastqc_data.txt"
        } else {
            copy_command = "cp ${prefix}_trimmed_paired_1_fastqc/summary.txt ${prefix}_1_trimmed_summary.txt && cp ${prefix}_trimmed_paired_2_fastqc/summary.txt ${prefix}_2_trimmed_summary.txt"
            data_command = "cp ${prefix}_trimmed_paired_1_fastqc/fastqc_data.txt ${prefix}_1_trimmed_fastqc_data.txt && cp ${prefix}_trimmed_paired_2_fastqc/fastqc_data.txt ${prefix}_2_trimmed_fastqc_data.txt"
        }
        '''
        !{params.fastqc} -t !{task.cpus} --extract !{params.fastqc_args} !{fastqc_trimmed_input}
        !{copy_command}
        !{data_command}
        '''
}

/*
 * Run the ERCC process if the --ercc flag is specified
 */

process ERCC {
	
    tag "$prefix"
    publishDir "${params.output}/ERCC/${prefix}",'mode':'copy'

    input:
        path ercc_index
        tuple val(prefix), path(ercc_input)
        path complete_manifest

    output:
        path "${prefix}_ercc_abundance.tsv", emit: ercc_abundances
        path "ercc_${prefix}.log"

    shell:
        if (params.sample == "single") {
            kallisto_flags = params.kallisto_quant_ercc_single_args
        } else {
            kallisto_flags = params.kallisto_quant_ercc_paired_args
        }
        '''
        ( set -o posix ; set ) > bash_vars.txt

        #  Find this sample's strandness and determine kallisto flags to set
        strand=$(cat samples_complete.manifest | grep " !{prefix} " | awk -F ' ' '{print $NF}')
        if [ $strand == 'forward' ]; then
            kallisto_strand=" --fr-stranded"
        elif [ $strand == 'reverse' ]; then
            kallisto_strand=" --rf-stranded"
        else
            kallisto_strand=""
        fi

        #  Quantify ERCC, and don't throw an error if 0 counts are found
        !{params.kallisto} quant \
            -i !{ercc_index} \
            -t !{task.cpus} \
            !{kallisto_flags} \
            -o . \
            $kallisto_strand \
            !{ercc_input} \
            || true

        cp abundance.tsv !{prefix}_ercc_abundance.tsv
        cp .command.log ercc_!{prefix}.log

        temp=$(( set -o posix ; set ) | diff bash_vars.txt - | grep '>' | cut -d ' ' -f 2- || true)
        echo "$temp" > bash_vars.txt
        '''
}
