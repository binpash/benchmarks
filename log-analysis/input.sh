#!/bin/bash

set -e

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/inputs"

mkdir -p "$input_dir"

for arg in "$@"; do
    if [[ "$arg" == "--min" ]]; then
        if [[ -d "$input_dir/nginx-logs-min" ]] && [[ -d "$input_dir/pcaps-min" ]]; then
            echo "Data already downloaded and extracted."
            exit 0
        fi
        cp -r "${eval_dir}/min_inputs/nginx-logs-min" "$input_dir"
        cp -r "${eval_dir}/min_inputs/pcaps-min" "$input_dir"
        exit 0
    fi
done

for arg in "$@"; do
    if [[ "$arg" == "--small" ]]; then
        url=https://atlas-group.cs.brown.edu/data/pcaps.zip
        zip_dst="$input_dir/pcaps.zip"
        wget --no-check-certificate $url -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        rm "$zip_dst"

        url=https://atlas-group.cs.brown.edu/data/nginx.zip
        zip_dst="$input_dir/nginx.zip"
        wget --no-check-certificate $url -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        rm "$zip_dst"

        exit 0
    fi
done

url=https://atlas-group.cs.brown.edu/data

zip_dst="$input_dir/pcaps.zip"
wget --no-check-certificate "$url"/pcaps.zip -O "$zip_dst"
unzip "$zip_dst" -d "$input_dir"
rm "$zip_dst"

wget --no-check-certificate "$url"/pcaps_large.zip -O "$zip_dst"
unzip "$zip_dst" -d "$input_dir"/pcaps
rm "$zip_dst"

zip_dst="$input_dir/nginx.zip"
wget --no-check-certificate "$url"/nginx.zip -O "$zip_dst"
unzip "$zip_dst" -d "$input_dir"
rm "$zip_dst"

# TODO: Fix the large nginx download
# wget --no-check-certificate "$url"/nginx_large.zip -O "$zip_dst"
# unzip "$zip_dst" -d "$input_dir/nginx"
# rm "$zip_dst"
