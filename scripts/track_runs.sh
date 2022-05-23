#   Script to log metadata about SPEAQeasy runs submitted at JHPCE. This helps
#   decentralize knowledge about which RNA-seq datasets have been processed,
#   where they were processed, and other info

log_dir='/dcs04/lieber/lcolladotor/libdData_LIBD001/processed-data/nextflow_logs/SPEAQeasy'
SPEAQeasy_log=$1

#   First, only continue if the SPEAQeasy log indicates completion of the
#   pipeline
if [[ $(grep -E "^Completed at: " $SPEAQeasy_log | wc -l) -ne 1 ]]; then
    exit 0
fi

#   Test runs should not be logged
is_small_test=$(
    grep -E "^Small test *: (false|true)$" $SPEAQeasy_log |
    tail -n 1 | cut -d ":" -f 2 | tr -d " "
)
if [[ "$is_small_test" == "true" ]]; then
    exit 0
fi

#   Grab the manifest path from the SPEAQeasy output log, then count the number
#   of samples
man_path=$(
    grep -E "^Input dir *: /.*" $SPEAQeasy_log |
    tail -n 1 | cut -d ":" -f 2 | tr -d " "
)/samples.manifest

pairedness=$(
    grep -E "^Sample *: (single|paired)$" $SPEAQeasy_log |
    tail -n 1 | cut -d ":" -f 2 | tr -d " "
)
if [[ "$pairedness" == "single" ]]; then
    num_samples=$(cut -f 3 $man_path | sort -u | wc -l)
elif [[ "$pairedness" == "paired" ]]; then
    num_samples=$(cut -f 5 $man_path | sort -u | wc -l)
else
    echo "While tracking run: could not determine if samples were paired-end."
    exit 1
fi

#   Tab-separated list:
#   [work dir] [date] [number of samples] [user]

log_path=$(mktemp -p $log_dir -t run_XXXXXXX.log)
echo -e "$PWD\t$(date +%Y-%m-%d,%H:%M)\t${num_samples}\t$(whoami)" > $log_path
chmod 755 $log_path
