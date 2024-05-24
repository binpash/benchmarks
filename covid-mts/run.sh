#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

if [[ "$1" == "--small" ]]; then
    echo "Using small input"
    input_file="$SUITE_DIR/inputs/in_small.csv"
else
    echo "Using default input"
    input_file="$SUITE_DIR/inputs/in.csv"
fi

mkdir -p "outputs"
all_res_file="./outputs/covid-mts.res"
> $all_res_file

# time_file stores the time taken for each script
# mode_res_file stores the time taken and the script name for every script in a mode (e.g. bash, pash, dish, fish)
# all_res_file stores the time taken for each script for every script run, making it easy to copy and paste into the spreadsheet
covid-mts() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/covid-mts.res"
    > $mode_res_file

    echo executing covid-mts $1 $(date) | tee -a $mode_res_file $all_res_file

    for number in `seq 4` ## initial: FIXME 5.sh is not working yet
    do
        script="${number}"
        script_file="./scripts/$script.sh"
        output_dir="./outputs/$1/$script/"
        output_file="./outputs/$1/$script.out"
        time_file="./outputs/$1/$script.time"
        log_file="./outputs/$1/$script.log"

        if [[ "$1" == "bash" ]]; then
            (time bash $script_file $input_file > $output_file ) 2> $time_file
        fi

        cat "${time_file}" >> $all_res_file
        echo "$script_file $(cat "$time_file")" | tee -a $mode_res_file
    done
}

covid-mts "bash"
