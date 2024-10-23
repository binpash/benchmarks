#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
input_dir="${eval_dir}/input"

mkdir -p $input_dir

data_url=https://atlas-group.cs.brown.edu/data/wav.zip
zip_dst=$input_dir/wav.zip
full_dir="$input_dir/wav_full"
small_dir="$input_dir/wav_small"
wget --no-check-certificate "$data_url" -O "$zip_dst"
unzip "$zip_dst" -d "$input_dir"
mkdir -p "$full_dir" "$small_dir"
# copy `.wav`s to their final destinations. 
# Make sure we have the correct number of inputs
# with numbered backups (do not overwrite inputs).
for i in {1..120}; do
    cp --backup=numbered $input_dir/wav/* "--target-directory=$full_dir"
done
for i in {1..10}; do
    cp --backup=numbered $input_dir/wav/* "--target-directory=$small_dir"
done
rm -r "$zip_dst" "$input_dir/wav"

data_url=https://atlas-group.cs.brown.edu/data/full/jpg.zip
zip_dst="$input_dir/jpg_full.zip"
out_dir="$input_dir/jpg_full"
wget --no-check-certificate $data_url -O $zip_dst
unzip $zip_dst -d $out_dir
rm "$zip_dst"

data_url=https://atlas-group.cs.brown.edu/data/small/jpg.zip
zip_dst=$input_dir/jpg_small.zip
out_dir=$input_dir/jpg_small
wget --no-check-certificate $data_url -O $zip_dst
unzip $zip_dst -d $out_dir
rm "$zip_dst"
