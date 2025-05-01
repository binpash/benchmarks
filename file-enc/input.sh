#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
input_dir="${eval_dir}/inputs"

mkdir -p "$input_dir"

DATA_LINK="https://atlas-group.cs.brown.edu/data"
ZIP_DST="$input_dir/pcaps.zip"

for arg in "$@"; do
    if [ "$arg" = "--min" ]; then
        mkdir -p "$input_dir/pcaps"
        cp min_inputs/* "$input_dir/pcaps/"
        exit 0
    fi
    if [ "$arg" = "--small" ]; then
        curl --insecure $DATA_LINK/pcaps.zip -o "$input_dir/pcaps.zip"
        unzip "$input_dir/pcaps.zip" -d "$input_dir"
        rm "$input_dir/pcaps.zip"
        exit 0
    fi
done

if [ ! -d "$input_dir/pcaps" ]; then
    mkdir -p "$input_dir"/pcaps

    wget --no-check-certificate "$DATA_LINK"/pcaps_large.zip -O "$ZIP_DST"
    unzip "$ZIP_DST" -d "$input_dir"/pcaps
    rm "$ZIP_DST"

    curl --insecure $DATA_LINK/pcaps.zip -o "$input_dir/pcaps.zip"
    unzip "$input_dir/pcaps.zip" -d "$input_dir"
    rm "$input_dir/pcaps.zip"
else
    echo "pcaps already exists. Skipping download."
fi
