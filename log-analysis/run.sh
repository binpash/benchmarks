#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

names_scripts=(
    "LogAnalysis1;nginx"
    "LogAnalysis2;pcaps"
  )

if [[ "$@" == *"--small"* ]]; then
    scripts_inputs=(
        "nginx;/log_data_small"
        "pcaps;/pcap_data_small"
    )
    scripts_outputs=(
        "nginx;/log-analysis/nginx_analysis_small"
        "pcaps;/log-analysis/pcap_analysis_small"
    )
else
    scripts_inputs=(
        "nginx;/log_data"
        "pcaps;/pcap_data"
    )
    scripts_outputs=(
        "nginx;/log-analysis/nginx_analysis"
        "pcaps;/log-analysis/pcap_analysis"
    )
fi

parse_directories() {
    local script_name=$1
    local scripts_array=("${!2}")
    for entry in "${scripts_array[@]}"; do
        IFS=";" read -r -a parsed <<< "${entry}"
        if [[ "${parsed[0]}" == "${script_name}" ]]; then
            echo "${parsed[1]}"
            return
        fi
    done
}

mkdir -p "outputs"
all_res_file="./outputs/log-analysis.res"
> $all_res_file

# time_file stores the time taken for each script
# mode_res_file stores the time taken and the script name for every script in a mode (e.g. bash, pash, dish, fish)
# all_res_file stores the time taken for each script for every script run, making it easy to copy and paste into the spreadsheet
log-analysis() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/log-analysis.res"
    > $mode_res_file

    echo executing log-analysis $1 $(date) | tee -a $mode_res_file $all_res_file
    for name_script in ${names_scripts[@]}
    do
        IFS=";" read -r -a name_script_parsed <<< "${name_script}"
        name="${name_script_parsed[0]}"
        script="${name_script_parsed[1]}"
        script_file="./scripts/$script.sh"
        input_dir=$(parse_directories "$script" scripts_inputs[@])
        output_dir=$(parse_directories "$script" scripts_outputs[@])
        output_file="./outputs/$1/$script.out"
        time_file="./outputs/$1/$script.time"
        log_file="./outputs/$1/$script.log"
        hash_file="./outputs/$1/$script.hash"


        if [[ "$1" == "bash" ]]; then
            (time $script_file $input_dir $output_dir > $output_file) 2> $time_file
        else
            params="$2"
            (time $PASH_TOP/pa.sh $params --log_file $log_file $script_file $input_dir $output_dir > $output_file) 2> $time_file
        fi

        # Generate SHA-256 hash and delete output file
        shasum -a 256 "$output_file" | awk '{ print $1 }' > "$hash_file"
        rm "$output_file"

        cat "${time_file}" >> $all_res_file
        echo "$script_file $(cat "$time_file")" | tee -a $mode_res_file
    done
}


# adjust the debug flag as required
d=0

log-analysis "bash"
log-analysis "pash"        "--width 8 --r_split -d $d --parallel_pipelines --profile_driven"
