#!/bin/bash

export TIMEFORMAT=%R
cd "$(realpath $(dirname "$0"))"
eval_dir="./scripts"

if [[ "$1" == "--small" ]]; then
    echo "Using small input"
    input_dir="./inputs/pcap_data_small"
else
    echo "Using default input"
    input_dir="./inputs/pcap_data"
fi


scripts=(
    "compress_files"
    "encrypt_files"
)

mkdir -p "outputs"
all_res_file="./outputs/file-enc.res"
> $all_res_file

file-enc() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/file-enc.res"
    > $mode_res_file

    echo executing file-enc $1 $(date) | tee -a $mode_res_file $all_res_file

    for script in ${scripts[@]}
    do
        IFS=";" read -r -a name_script_parsed <<< "${name_script}"
        script_file="./scripts/$script.sh"
        # input for all nlp scripts is ./inputs/pg, which is already default for each script
        output_dir="./outputs/$1/$script/"
        output_file="./outputs/$1/$script.out"
        time_file="./outputs/$1/$script.time"
        log_file="./outputs/$1/$script.log"

        # output_file contains "done" when run successfully. The real outputs are under output_dir/
        if [[ "$1" == "bash" ]]; then
            (time bash $script_file $input_dir $output_dir > $output_file ) 2> $time_file
        else
            params="$2"
            if [[ $2 == *"--ft"* ]]; then
                params="$2 --script_name $script_file"
            fi

            (time $PASH_TOP/pa.sh $params --log_file $log_file $script_file $input_dir $output_dir > $output_file) 2> $time_file

            if [[ $2 == *"--kill"* ]]; then
                sleep 10
                python3 "$DISH_TOP/evaluation/notify_worker.py" resurrect
            fi

            sleep 10
        fi

        cat "${time_file}" >> $all_res_file
        echo "$script_file $(cat "$time_file")" | tee -a $mode_res_file
    done
}

d=0

file-enc "bash"
file-enc "pash"        "--width 8 --r_split --parallel_pipelines --profile_driven -d $d"