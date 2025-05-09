#!/bin/bash

export PATH=$PATH:$HOME/.local/bin
export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
export BIODIR="${PASH_SPEC_TOP}/report/benchmarks/bio4"

download_dir="${PASH_SPEC_TOP}/report/resources/bio4/"
echo $download_dir
mkdir -p "$download_dir"
INPUT_LIST=$1
input_file="${BIODIR}/$INPUT_LIST"

sudo apt install -y samtools time

download_bam_files() {
    local input_file=$1
    while IFS='-' read -r pop sample link; do
        if [[ ! -f "$download_dir/$sample.bam" ]]; then
            echo "Downloading $sample..."
            wget -q -O "$download_dir/$sample.bam" "$link"
        fi
    done < "$input_file"
}

while getopts ":i:c" opt; do
    case $opt in
        i) input_file=$OPTARG ;;
        c) 
            echo "Cleaning up..."
            rm -rf "$download_dir"/*.bam
            rm -rf "$download_dir"/*.sam
            rm -rf "$download_dir"/../output
            exit 0
            ;;
        \?) echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
    esac
done

download_bam_files "$input_file"
