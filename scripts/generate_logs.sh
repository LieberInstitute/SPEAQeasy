#!/bin/bash

#  The full path to the main output log from nextflow is passed to this script
LOG=$1

#  Infer working directory and sample names from the main output log
work_dir=$(grep "^Working dir" $LOG | head -n 1 | tr -d [:blank:] | cut -d ":" -f 2)
man=$(grep "^Input" $LOG | head -n 1 | tr -d [:blank:] | cut -d ":" -f 2)/samples.manifest
samp_names=$(awk '{print $NF}' $man | uniq)

#  Place command logs in the output directory
out_dir=$(grep "^Output dir" $LOG | head -n 1 | tr -d [:blank:] | cut -d ":" -f 2)
mkdir -p $out_dir/logs

for samp in $samp_names; do
    i=1
    out_file="$out_dir/logs/${samp}_process_trace.log"
    rm -f $out_file
    echo "#########################" >> $out_file
    echo "  Sample $samp" >> $out_file
    echo -e "#########################\n\n" >> $out_file
    
    #  Grab all unique process submissions
    proc_lines=$(grep "\[.*/.*\] process.*\>" $LOG | cut -d "[" -f 1-2 | awk '!seen[$0]++')
    
    #  Loop through these
    for line_num in `seq 1 $(echo "$proc_lines" | wc -l)`; do
        proc_line=$(echo "$proc_lines" | awk "NR==$line_num")
        
        #  Infer process working directory
        suffix=$(echo $proc_line | cut -d '[' -f 2 | cut -d ']' -f 1)
        proc_dir=$(echo $work_dir/${suffix}*)
        
        #  If files from this sample can be found in the working directory
        if [ $(ls $proc_dir | grep $samp | wc -l) -gt 0 ]; then
        
            #  Infer process name and exit code
            proc_name=$(echo $proc_line | cut -d " " -f 4)
            exit_code=$(cat $proc_dir/.exitcode)
            
            #  Print the process name, working directory, exit code, and command run
            echo "######################################################" >> $out_file
            echo "  Process $i: ${proc_name}" >> $out_file
            echo -e "######################################################\n" >> $out_file
            echo "    Working directory: ${proc_dir}" >> $out_file
            echo "    Exit code for process: ${exit_code}" >> $out_file         
            echo -e "    Command run:\n" >> $out_file
            echo "--------------------------------- BEGIN COMMANDS -------------" >> $out_file
            cd $proc_dir
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
            
            i=$((++i))
        fi
    done
    
    echo "This was the last process submitted for sample $samp." >> $out_file
done