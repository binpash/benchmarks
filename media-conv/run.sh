#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
export TIMEFORMAT=%R
cd $SUITE_DIR

names_scripts=(
    "MediaConv1;img_convert"
    "MediaConv2;to_mp3"
  )

if [[ "$@" == *"--small"* ]]; then
    scripts_inputs=(
        "img_convert;/jpg_small/jpg"
        "to_mp3;/wav_small/wav"
    )
    scripts_outputs=(
        "img_convert;jpg_small"
        "to_mp3;mp3_small"
    )
else
    scripts_inputs=(
        "img_convert;/jpg_full/jpg"
        "to_mp3;/wav_full/wav"
    )
    scripts_outputs=(
        "img_convert;jpg"
        "to_mp3;mp3"
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
all_res_file="./outputs/media-conv.res"
> $all_res_file

gen-hash-dir() {
    # Iterate over each file in the directory and subdirectories, excluding any .hash files
    find "$1" -type f ! -name "*.hash" | while read -r file; do
    if [ -f "$file" ]; then
        # Compute the SHA-256 hash of the file
        hash_value=$(sha256sum "$file" | awk '{print $1}')
        
        # Extract the basename without the extension
        base_name=$(basename "$file" .zip)
        
        # Create the .hash filename
        hash_file="${1}/${base_name}.hash"
        ls 
        # Write the hash to the .hash file
        echo "$hash_value" > "$hash_file"
        
        # Print a message indicating the hash file has been created
        # echo "Hash for $(basename "$file") written to ${base_name}.hash"
    fi
    done
}

# time_file stores the time taken for each script
# mode_res_file stores the time taken and the script name for every script in a mode (e.g. bash, pash, dish, fish)
# all_res_file stores the time taken for each script for every script run, making it easy to copy and paste into the spreadsheet
media-conv() {
    mkdir -p "outputs/$1"
    mode_res_file="./outputs/$1/media-conv.res"
    > $mode_res_file

    echo executing media-conv $1 $(date) | tee -a $mode_res_file $all_res_file
    for name_script in ${names_scripts[@]}
    do
        IFS=";" read -r -a name_script_parsed <<< "${name_script}"
        name="${name_script_parsed[0]}"
        script="${name_script_parsed[1]}"
        script_file="./scripts/$script.sh"
        input_dir=$(parse_directories "$script" scripts_inputs[@])
        output_dir=./outputs/$1/$(parse_directories "$script" scripts_outputs[@])
        output_file="./outputs/$1/$script.out"
        time_file="./outputs/$1/$script.time"
        log_file="./outputs/$1/$script.log"
        hash_file="./outputs/$1/$script.hash"


        if [[ "$1" == "bash" ]]; then
            (time $script_file $input_dir $output_dir > $output_file) 2> $time_file
        fi

        # Generate hash files for each converted file in output_dir
        gen-hash-dir "$output_dir"
        
        # Find and delete all files not ending with 'hash'
        find "$output_dir" -type f ! -name '*hash' -delete

        cat "${time_file}" >> $all_res_file
        echo "$script_file $(cat "$time_file")" | tee -a $mode_res_file
    done
}


# adjust the debug flag as required
d=0

media-conv "bash"

