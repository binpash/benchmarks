#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
input_dir="${eval_dir}/input"
mkdir -p "$input_dir"

# dpt source images
data_url="https://atlas-group.cs.brown.edu/data/dpt.zip"
zip_dst="${input_dir}/images.zip"
wget --no-check-certificate "$data_url" -O "$zip_dst"
unzip "$zip_dst" -d "$input_dir"
rm "$zip_dst"

full_dir="${input_dir}/images_full"
small_dir="${input_dir}/images_small"
mkdir -p "$full_dir" "$small_dir"

for i in {1..120}; do
    cp --backup=numbered "$input_dir"/images/* "$full_dir"
done

for i in {1..10}; do
    cp --backup=numbered "$input_dir"/images/* "$small_dir"
done

rm -r "$input_dir/images"
