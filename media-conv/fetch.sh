#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
URL='https://atlas.cs.brown.edu/data'
input_dir="${REPO_TOP}/media-conv/inputs"

mkdir -p "$input_dir"

data_url=$URL/wav.zip
zip_dst=$input_dir/wav.zip
full_dir="$input_dir/wav_full"
small_dir="$input_dir/wav_small"
min_dir="$input_dir/wav_min"

if [[ ! -d "$full_dir" ]]; then
    wget --no-check-certificate "$data_url" -O "$zip_dst"
    unzip "$zip_dst" -d "$input_dir"
    rm -rf "$full_dir" "$small_dir" "$min_dir"
    mkdir -p "$full_dir" "$small_dir" "$min_dir"

    # copy `.wav`s to their final destinations.
    # Make sure we have the correct number of inputs
    # with numbered backups (do not overwrite inputs).
    for i in {1..120}; do
        cp --backup=numbered "$input_dir"/wav/* "--target-directory=$full_dir"
    done
    for i in {1..10}; do
        cp --backup=numbered "$input_dir"/wav/* "--target-directory=$small_dir"
    done
    for i in {1..2}; do
        cp --backup=numbered "$input_dir"/wav/* "--target-directory=$min_dir"
    done
    rm -r "$zip_dst" "$input_dir/wav"
fi

# if small flag
for arg in "$@"; do
    if [[ "$arg" == "--small" ]]; then
        if [[ -d "$input_dir/jpg_small" ]]; then
            echo "Data already downloaded and extracted."
            exit 0
        fi
        data_url=${URL}/small/jpg.zip
        zip_dst=$input_dir/jpg_small.zip
        out_dir=$input_dir/jpg_small
        wget --no-check-certificate $data_url -O $zip_dst || {
            echo "Failed to download $data_url"
            exit 1
        }
        unzip $zip_dst -d $out_dir || {
            echo "Failed to unzip $zip_dst"
            exit 1
        }
        rm "$zip_dst"
        exit 0
    elif [[ "$arg" == "--min" ]]; then
        if [[ -d "$input_dir/jpg_min" ]]; then
            echo "Data already downloaded and extracted."
            exit 0
        fi
        min_inputs="${REPO_TOP}/media-conv/min_inputs/"
        mkdir -p "$input_dir"
        cp -r "$min_inputs"/* "$input_dir/"
        exit 0
    fi
done

if [[ -d "$input_dir/jpg" ]]; then
    echo "Data already downloaded and extracted."
    exit 0
fi
data_url=${URL}/full/jpg.zip
zip_dst="$input_dir/jpg_full.zip"
out_dir="$input_dir/jpg_full"
wget --no-check-certificate $data_url -O $zip_dst
unzip $zip_dst -d $out_dir
rm "$zip_dst"
