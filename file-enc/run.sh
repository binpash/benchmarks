#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
input_dir="${eval_dir}/input"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
mkdir -p $results_dir

input_pcaps="$input_dir/pcaps"
suffix=".full"
if [[ "$1" == "--small" ]]; then
    # TODO: prepare a smaller input
    input_pcaps="$input_dir/pcaps"
    suffix=".small"
fi

$scripts_dir/compress_files.sh $input_pcaps $results_dir/compress_files$suffix
$scripts_dir/encrypt_files.sh $input_pcaps $results_dir/encrypt_files$suffix
