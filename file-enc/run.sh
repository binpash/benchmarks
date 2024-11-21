#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
input_dir="${eval_dir}/input"
results_dir="${eval_dir}/results"
scripts_dir="${eval_dir}/scripts"
mkdir -p $results_dir

echo executing file-enc $(date)

input_pcaps="$input_dir/pcaps"
suffix=".full"
if [[ "$1" == "--small" ]]; then
    # TODO: prepare a smaller input
    input_pcaps="$input_dir/pcaps"
    suffix=".small"
fi

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
$BENCHMARK_SHELL $scripts_dir/compress_files.sh $input_pcaps $results_dir/compress_files$suffix
$BENCHMARK_SHELL $scripts_dir/encrypt_files.sh $input_pcaps $results_dir/encrypt_files$suffix
