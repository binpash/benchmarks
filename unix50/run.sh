#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

if [[ "$@" == *"--small"* ]]; then
    scripts_inputs=(
        "1;1_1M"
        "2;1_1M"
        "3;1_1M"
        "4;1_1M"
        "5;2_1M"
        "6;3_1M"
        "7;4_1M"
        "8;4_1M"
        "9;4_1M"
        "10;4_1M"
        "11;4_1M"
        "12;4_1M"
        "13;5_1M"
        "14;6_1M"
        "15;7_1M"
        "16;7_1M"
        "17;7_1M"
        "18;8_1M"
        "19;8_1M"
        "20;8_1M"
        "21;8_1M"
        # "22;8_1M"
        "23;9.1_1M"
        "24;9.2_1M"
        "25;9.3_1M"
        "26;9.4_1M"
        # "27;9.5_1M"
        "28;9.6_1M"
        "29;9.7_1M"
        "30;9.8_1M"
        "31;9.9_1M"
        "32;10_1M"
        "33;10_1M"
        "34;10_1M"
        "35;11_1M"
        "36;11_1M"
    )
else
        scripts_inputs=(
        "1;1_3G"
        "2;1_3G"
        "3;1_3G"
        "4;1_3G"
        "5;2_3G"
        "6;3_3G"
        "7;4_3G"
        "8;4_3G"
        "9;4_3G"
        "10;4_3G"
        "11;4_3G"
        "12;4_3G"
        "13;5_3G"
        "14;6_3G"
        "15;7_3G"
        "16;7_3G"
        "17;7_3G"
        "18;8_3G"
        "19;8_3G"
        "20;8_3G"
        "21;8_3G"
        # "22;8_3G"
        "23;9.1_3G"
        "24;9.2_3G"
        "25;9.3_3G"
        "26;9.4_3G"
        # "27;9.5_3G"
        "28;9.6_3G"
        "29;9.7_3G"
        "30;9.8_3G"
        "31;9.9_3G"
        "32;10_3G"
        "33;10_3G"
        "34;10_3G"
        "35;11_3G"
        "36;11_3G"
    )
fi

mkdir -p "outputs"
all_res_file="./outputs/unix50.res"
> $all_res_file

# time_file stores the time taken for each script
# mode_res_file stores the time taken and the script name for every script in a mode (e.g. bash, pash, dish, fish)
# all_res_file stores the time taken for each script for every script run, making it easy to copy and paste into the spreadsheet
unix50() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/unix50.res"
    > $mode_res_file

    echo executing unix50 $1 $(date) | tee -a $mode_res_file $all_res_file

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

unix50 "bash"
