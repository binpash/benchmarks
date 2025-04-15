#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/image_annotation"
inputs_dir="$eval_dir/inputs"

mkdir -p "$inputs_dir"
flag=0
suffix=""
size=100000000
for arg in "$@"; do
    if [[ "$arg" == "--small" ]]; then
        data_url=https://atlas-group.cs.brown.edu/data/small/jpg.zip
        zip_dst=$inputs_dir/jpg.small.zip
        out_dir=$inputs_dir/jpg.small
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
        min_inputs="$eval_dir/min_inputs/"
        out_dir="$inputs_dir/jpg.min"
        mkdir -p "$out_dir"
        cp -r "$min_inputs"/* "$out_dir/"
        exit 0
    fi
done

echo "Downloading full dataset."
data_url=https://atlas-group.cs.brown.edu/data/full/jpg.zip
zip_dst="$inputs_dir/jpg.zip"
out_dir="$inputs_dir/jpg"
wget --no-check-certificate $data_url -O $zip_dst
unzip $zip_dst -d $out_dir
rm "$zip_dst"
