#!/bin/bash

set -e

# creates input/pcaps and input/nginx-logs

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/inputs"

mkdir -p "$input_dir"

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
