#   Script to log metadata about SPEAQeasy runs submitted at JHPCE. This helps
#   decentralize knowledge about which RNA-seq datasets have been processed,
#   where they were processed, and other info

log_dir='/dcs04/lieber/lcolladotor/libdData_LIBD001/processed-data/nextflow_logs/SPEAQeasy'
SPEAQeasy_log='SPEAQeasy_output.log'

#   First, only continue if the SPEAQeasy log indicates completion of the
#   pipeline
if [[ $(grep -E "^Completed at: " $SPEAQeasy_log | wc -l) -ne 1 ]]; then
    exit 0
fi

#   Test runs should not be logged
is_small_test=$(grep -E "^Small test *: (false|true)$" $SPEAQeasy_log | cut -d ":" -f 2 | tr -d " ")
if [[ "$is_small_test" == "true" ]]; then
    exit 0
fi

#   Grab the manifest path from the SPEAQeasy output log, then count the number
#   of samples
man_path=$(grep -E "^Input dir *: /.*" $SPEAQeasy_log | cut -d ":" -f 2 | tr -d " ")/samples.manifest
num_samples=$(cut -f 5 $man_path | sort -u | wc -l)

#   Tab-separated list:
#   [work dir] [date] [number of samples] [user]

log_path=$(mktemp -p $log_dir)
echo -e "$PWD\t$(date +%Y-%m-%d,%H:%M)\t${num_samples}\t$(whoami)" > $log_path
