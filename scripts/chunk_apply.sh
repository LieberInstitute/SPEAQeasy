#   This script is intended to a apply a given command in chunks. In particular,
#   it can take any command that works with N files, and it applies it in chunks
#   of B files, iteratively, until there are at most B files left. The
#   motivation is that there are several tasks in SPEAQeasy (computing coverage,
#   merging variant calls) that can fail with sufficiently many samples. Reasons
#   for failure include system limits on number of open file handles or threads.
#   This script can break arbitrarily large datasets into sufficiently small
#   batches to avoid these errors. 
#
#   Tested this script with a working directory of different numbers of
#   'temp.[unique number].txt' files, different batch sizes, and:
#       file_regex='.*\.txt$'
#       command="cat [files] > temp_auto_[index].txt"

#  Rename command-line arguments
batch_size=$1
file_regex=$2
command=$3

file_list=$(ls -1 | grep -E "$file_regex")
num_files=$(echo "$file_list" | wc -l)
num_batches=$(($num_files / $batch_size))
iter_num=0

while [ $num_files -gt $batch_size ]; do
    echo "Dividing existing files into batches of size ${batch_size}, which involves ${num_batches} batches..."
    for i in $(seq 1 $num_batches); do
        #  Form the batch of files
        start_index=$((($i - 1) * $batch_size + 1))
        end_index=$(($i * $batch_size))
        temp_files=$(echo "$file_list" | sed -n "${start_index},${end_index}p" | paste -sd " ")

        #  Apply the command to the batch and output a unique temporary file.
        #  Note that the use of 'eval' here is fairly safe, and can only have
        #  harmful results if params.bcftools, params.wiggletools, or
        #  params.bgzip are intentionally set to something crazy
        eval $(echo $command | sed "s%\[index\]%${iter_num}_$i%" | sed "s%\[files\]%${temp_files}%")
        rm $temp_files
    done
    
    #  If the remaining piece exists and involves more than one file
    if [ $num_files -gt $(($num_batches * $batch_size + 1)) ]; then
        echo "On remaining piece..."
        start_index=$(($num_batches * $batch_size + 1))
        temp_files=$(echo "$file_list" | sed -n "${start_index},${num_files}p" | paste -sd " ")
        
        #  Apply the command to the remaining piece and replace it with a unique
        #  temporary file
        eval $(echo $command | sed "s%\[index\]%${iter_num}_$(($num_batches + 1))%" | sed "s%\[files\]%${temp_files}%")
        rm $temp_files
    fi
    
    file_list=$(ls -1 | grep -E "$file_regex")
    num_files=$(echo "$file_list" | wc -l)
    num_batches=$(($num_files / $batch_size))
    iter_num=$(($iter_num + 1))
done
