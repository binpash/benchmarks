#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd "$SUITE_DIR"

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
        "shortest-scripts;all_cmds"
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
        "shortest-scripts;all_cmdsx100"
    )
fi

mkdir -p "outputs"

oneliners() {
    echo executing oneliners $1 $(date)

    mkdir -p "outputs/$1"

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

        echo "$script_file $(cat "$time_file")" 
    done
}

oneliners "bash"
