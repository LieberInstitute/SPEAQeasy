#  Rename command-line arguments
sample=$1
num_reads_infer_strand=$2
kallisto=$3
num_threads=$4
kallisto_len_mean=$5
kallisto_len_sd=$6
strand=$7
rscript=$8
strand_mode=$9
id=${10}

if [ $num_threads == 1 ]; then
    thread_opt=""
else
    thread_opt="-t $num_threads"
fi
        
if [ $sample == "paired" ]; then
    fq1=$(ls *_1.f*q*)
    fq2=$(ls *_2.f*q*)

    #  Subset the FASTQ files by the given number of reads
    if [ $(basename $fq1 | grep ".gz" | wc -l) == 1 ]; then
        zcat $fq1 | head -n $num_reads_infer_strand > test_1.fastq
        zcat $fq2 | head -n $num_reads_infer_strand > test_2.fastq
    else
        cat $fq1 | head -n $num_reads_infer_strand > test_1.fastq
        cat $fq2 | head -n $num_reads_infer_strand > test_2.fastq
    fi
        
    #  Try pseudoalignment assuming each type of strandness possible
    echo "Testing pseudoalignment assuming forward-strandness..."
    $kallisto quant $thread_opt --fr-stranded -i kallisto_index_* -o ./forward test_1.fastq test_2.fastq
    echo "Testing pseudoalignment assuming reverse-strandness..."
    $kallisto quant $thread_opt --rf-stranded -i kallisto_index_* -o ./reverse test_1.fastq test_2.fastq
else
    fq=$(ls *.f*q*)
            
    #  Subset the FASTQ file by the given number of reads
    if [ $(basename $fq | grep ".gz" | wc -l) == 1 ]; then
        zcat $fq | head -n $num_reads_infer_strand > test.fastq
    else
        cat $fq | head -n $num_reads_infer_strand > test.fastq
    fi
            
    #  Try pseudoalignment assuming each type of strandness possible
    echo "Testing pseudoalignment assuming forward-strandness..."
    $kallisto quant \
        $thread_opt \
        --single \
        -l $kallisto_len_mean \
        -s $kallisto_len_sd \
        --fr-stranded \
        -i kallisto_index_* \
        -o ./forward test.fastq
    echo "Testing pseudoalignment assuming reverse-strandness..."
    $kallisto quant \
        $thread_opt \
        --single \
        -l $kallisto_len_mean \
        -s $kallisto_len_sd \
        --rf-stranded \
        -i kallisto_index_* \
        -o ./reverse test.fastq
fi
        
forward_count=$(grep "n_pseudoaligned" forward/run_info.json | cut -d ":" -f 2 | cut -d "," -f 1)
reverse_count=$(grep "n_pseudoaligned" reverse/run_info.json | cut -d ":" -f 2 | cut -d "," -f 1)
        
#  Pass number of reads pseudoaligned for each test to the R script, which will
#  ultimately infer strandness (and print additional info/ potential warnings)
        
echo "Passing this info to the R script to infer strand..."
$rscript infer_strand.R -f $forward_count -r $reverse_count -p $sample -s $strand -m $strand_mode -i $id
if [ ! "$?" == 0 ]; then
    exit 1
else
    echo "Done."
fi
