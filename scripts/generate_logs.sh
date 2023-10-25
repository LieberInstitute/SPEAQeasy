#!/bin/bash

#  The full path to the main output log from nextflow is passed to this script
LOG=$1

#  Infer working directory and sample names from the main output log
work_dir=$(grep "^Working dir" $LOG | head -n 1 | tr -d [:blank:] | cut -d ":" -f 2)

in_dir=$(grep "^Input" $LOG | head -n 1 | tr -d [:blank:] | cut -d ":" -f 2)
if [[ -z $in_dir ]]; then
    echo "Incomplete SPEAQeasy log: not generating sample-specific logs."
    exit 1
fi
man=$in_dir/samples.manifest
samp_names=$(awk '{print $NF}' $man | uniq)

#  Place command logs in the output directory
out_dir=$(grep "^Output dir" $LOG | head -n 1 | tr -d [:blank:] | cut -d ":" -f 2)
mkdir -p $out_dir/logs

#  For pipeline runs halting with an error, determine if the failing process
#  is missing from the list of process submissions. For successful pipeline
#  runs, return false.
is_missing () {
    if [ $(grep -n "Work dir:" $LOG | wc -l) -gt 0 ]; then
        local line_num=$(($(grep -n "Work dir:" $LOG | head -n 1 | cut -d ":" -f 1) + 1))
        local suffix=$(sed -n "${line_num}p" $LOG | awk -F "/" '{print $NF}' | cut -c 1-6)
        if [ $(grep "$suffix]" $LOG | wc -l) -gt 0 ]; then
            return 0
        else
            return 1
        fi
    else
        return 0
    fi
}

print_commands () {
    local proc_dir=$1
    local proc_name=$2
    local exit_code=$3
    local samp=$4
    local out_file=$5
    local temp_i=$6
    
    #  If files from this sample can be found in the working directory
    if [ $(ls $proc_dir | grep $samp | wc -l) -gt 0 ]; then
        #  Print the process name, working directory, exit code, and command run
        echo "######################################################" >> $out_file
        echo "  Process $i: ${proc_name}" >> $out_file
        echo -e "######################################################\n" >> $out_file
        echo "    Working directory: ${proc_dir}" >> $out_file
        echo "    Exit code for process: ${exit_code}" >> $out_file         
        echo -e "    Command run:\n" >> $out_file
        echo "--------------------------------- BEGIN COMMANDS -------------" >> $out_file
        cd $proc_dir
        if [[ -f .command.sh ]]; then
          for j in `seq 1 $(cat .command.sh | wc -l)`; do
            #  Grab one line of the commands
            temp_line=$(awk "NR == $j" .command.sh)
                
            #  Loop through all defined bash variables, and substitute in the
            #  value of each variable appearing in the form "$variable_name".
            #  The second condition is a heuristic, intended to in practice
            #  determine if the process completed (bash_vars.txt has unintended
            #  content for incomplete processes, and the file tends to have
            #  many lines)
            if [ -f bash_vars.txt ] && [ 25 -gt $(cat bash_vars.txt | wc -l) ]; then
                for k in `seq 1 $(wc -l bash_vars.txt | cut -d " " -f 1)`; do
                    var_name=$(awk "NR == $k" bash_vars.txt | cut -d '=' -f 1)
                    var_val=$(awk "NR == $k" bash_vars.txt | cut -d '=' -f 2 | tr -d "'")
                        
                    temp_line=$(echo $temp_line | sed "s@\$${var_name}\b@${var_val}@g" | sed "s@\${$var_name}@$var_val@g")
                done
            fi
                
            #  The result is what is printed for the user in the log
            echo $temp_line >> $out_file
          done
        fi
        echo -e "--------------------------------- END COMMANDS ---------------\n" >> $out_file
            
        #  Link to log if process had exit code 0, otherwise warn and print log output
        if [ "$exit_code" == "0" ]; then
            echo -e "    Process appeared to complete without errors. Log is here: $proc_dir/.command.log\n\n" >> $out_file
        else
            echo "    Process had nonzero exit status (is it still running? Otherwise something may" >> $out_file
            echo -e "    have went wrong). The full log output is here:\n" >> $out_file
            echo "--------------------------------- BEGIN LOG -------------" >> $out_file
            cat $proc_dir/.command.log >> $out_file
            echo -e "--------------------------------- END LOG ---------------\n\n" >> $out_file
        fi
        temp_i=$(($temp_i + 1))
    fi
    
    return $temp_i
}

for samp in $samp_names; do
    i=1
    out_file="$out_dir/logs/${samp}_process_trace.log"
    rm -f $out_file
    echo "#########################" >> $out_file
    echo "  Sample $samp" >> $out_file
    echo -e "#########################\n\n" >> $out_file
    
    #  Grab all unique process submissions
    proc_lines=$(grep "\[.*/.*\] process.*\>" $LOG | cut -d "[" -f 1-2 | awk '!seen[$0]++')
    
    #  Loop through these, if there are any
    if [[ ! -z $proc_lines ]]; then
        total_processes=$(($(echo "$proc_lines" | wc -l) + 1))
        for line_num in $(seq 1 $total_processes); do
            #  The last process submission to check is an exceptional case: determined
            #  by the error message at the end of the log (which strangely sometimes
            #  lacks a typical process submission line above)
            if [ $line_num -eq $total_processes ]; then
                is_missing
                if [ $? -eq 1 ]; then
                    line_num_temp=$(($(grep -n "Work dir:" $LOG | head -n 1 | cut -d ":" -f 1) + 1))
                    
                    #  Infer process working directory, name, and exit code
                    proc_dir=$(sed -n "${line_num_temp}p" $LOG | tr -d " ")
                    proc_name=$(grep "Error executing process" $LOG | head -n 1 | cut -d '>' -f 2 | tr -d "'| ")
                    if [[ -f $proc_dir/.exitcode ]]; then
                      exit_code=$(cat $proc_dir/.exitcode)
                      print_commands $proc_dir $proc_name $exit_code $samp $out_file $i
                      i=$(echo $?)
                    fi
                fi
            else
                proc_line=$(echo "$proc_lines" | awk "NR==$line_num")
                suffix=$(echo $proc_line | cut -d '[' -f 2 | cut -d ']' -f 1)
                
                #  Infer process working directory, name, and exit code
                proc_dir=$(echo $work_dir/${suffix}*)
                proc_name=$(echo $proc_line | cut -d " " -f 4)

                if [[ -f $proc_dir/.exitcode ]]; then
                    exit_code=$(cat $proc_dir/.exitcode)
                else
                    exit_code=NA
                fi
                
                print_commands $proc_dir $proc_name $exit_code $samp $out_file $i
                i=$(echo $?)
            fi
        done
        
        echo "This was the last process submitted for sample $samp." >> $out_file
    fi
done
