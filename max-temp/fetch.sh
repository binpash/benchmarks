#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)

eval_dir="${REPO_TOP}/max-temp"
input_dir="${eval_dir}/inputs"

URL='https://atlas.cs.brown.edu/data'
URL=$URL/max-temp
FROM=2000
TO=2015

n_samples=99999
suffix="full"

mkdir -p "${input_dir}"

for arg in "$@"; do
    if [[ "$arg" == "--min" ]]; then
        if [[ -f "$input_dir/temperatures.min.txt" ]]; then
            echo "Data already downloaded and extracted."
            exit 0
        fi
        min_inputs="$eval_dir/min_inputs/"
        mkdir -p "$input_dir"
        cp -r "$min_inputs"/* "$input_dir/"
        exit 0
    fi
    if [[ "$arg" == "--small" ]]; then
        if [[ -f "$input_dir/temperatures.small.txt" ]]; then
            echo "Data already downloaded and extracted."
            exit 0
        fi
        data_url="${URL}/temperatures.small.tar.gz"
        wget --no-check-certificate "$data_url" -O "$input_dir/temperatures.small.tar.gz" || {
            echo "Failed to download $data_url"
            exit 1
        }
        tar -xzf "$input_dir/temperatures.small.tar.gz" -C "$input_dir" || {
            echo "Failed to extract $input_dir/temperatures.small.tar.gz"
            exit 1
        }
        rm "$input_dir/temperatures.small.tar.gz"
        mv "$input_dir/inputs/temperatures.small.txt" "$input_dir/temperatures.small.txt"
        rm -rf "$input_dir/inputs"
        exit 0
    fi
done
if [[ -f "$input_dir/temperatures.full.txt" ]]; then
    echo "Data already downloaded and extracted."
    exit 0
fi
data_url="${URL}/temperatures.full.tar.gz"
wget --no-check-certificate "$data_url" -O "$input_dir/temperatures.full.tar.gz" || {
    echo "Failed to download $data_url"
    exit 1
}
tar -xzf "$input_dir/temperatures.full.tar.gz" -C "$input_dir" || {
    echo "Failed to extract $input_dir/temperatures.full.tar.gz"
    exit 1
}
rm "$input_dir/temperatures.full.tar.gz"
mv "$input_dir/inputs/temperatures.full.txt" "$input_dir/temperatures.full.txt"
rm -rf "$input_dir/inputs"
exit 0
