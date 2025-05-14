#!/bin/bash
KOALA_SHELL=${KOALA_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/llm/scripts/image-annotation"
inputs_dir="$eval_dir/inputs"

mkdir -p "$inputs_dir"

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done


if [[ "$size" == "small" ]]; then
    # if inputs exist
    if [[ -d "$inputs_dir/jpg.small" ]]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi
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
elif [[ "$size" == "min" ]]; then
    if [[ -d "$inputs_dir/jpg.min" ]]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi
    min_inputs="$eval_dir/min_inputs/"
    out_dir="$inputs_dir/jpg.min/jpg"
    mkdir -p "$out_dir"
    cp -r "$min_inputs"/* "$out_dir/"
    exit 0
else
    if [[ -d "$inputs_dir/jpg" ]]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi

    echo "Downloading full dataset."
    data_url=https://atlas-group.cs.brown.edu/data/full/jpg.zip
    zip_dst="$inputs_dir/jpg.zip"
    out_dir="$inputs_dir/jpg"
    wget --no-check-certificate $data_url -O $zip_dst
    unzip $zip_dst -d $out_dir
    rm "$zip_dst"
fi