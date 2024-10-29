#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

if [[ "$@" == *"--small"* ]]; then
    scripts_inputs=(
        "nfa-regex;1M"
        "sort;1M"
        "top-n;1M"
        "wf;1M"
        "spell;1M"
        "diff;1M"
        "bi-grams;1M"
        "set-diff;1M"
        "sort-sort;1M"
    )
else
    scripts_inputs=(
        "nfa-regex;1G"
        "sort;3G"
        "top-n;3G"
        "wf;3G"
        "spell;3G"
        "diff;3G"
        "bi-grams;3G"
        "set-diff;3G"
        "sort-sort;3G"
    )
fi

mkdir -p "outputs"
all_res_file="./outputs/oneliners.res"
> $all_res_file

# time_file stores the time taken for each script
# mode_res_file stores the time taken and the script name for every script in a mode (e.g. bash, pash, dish, fish)
# all_res_file stores the time taken for each script for every script run, making it easy to copy and paste into the spreadsheet
oneliners() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/oneliners.res"
    > $mode_res_file

    echo executing oneliners $1 $(date) | tee -a $mode_res_file $all_res_file

    for script_input in ${scripts_inputs[@]}
    do
        IFS=";" read -r -a parsed <<< "${script_input}"
        script_file="./scripts/${parsed[0]}.sh"
        input_file="./inputs/${parsed[1]}.txt"
        output_file="./outputs/$1/${parsed[0]}.out"
        time_file="./outputs/$1/${parsed[0]}.time"
        log_file="./outputs/$1/${parsed[0]}.log"

        if [[ "$1" == "bash" ]]; then
            (time $script_file $input_file > $output_file) 2> $time_file
        fi

        cat "${time_file}" >> $all_res_file
        echo "$script_file $(cat "$time_file")" | tee -a $mode_res_file
    done
}

oneliners "bash"
